//
// Created by ylf9811 on 2021/7/6.
//
#include <functional>
#include <cassert>

#include "seqc.h"

SeQc::SeQc(CmdInfo *cmd_info1) {
    cmd_info = cmd_info1;
}

SeQc::~SeQc() {}

void SeQc::StateInfo(neoReference &ref) {
    int slen = ref.lseq;
    int qlen = ref.lqual;
    assert(slen == qlen);
    char *bases = reinterpret_cast<char *>(ref.base + ref.pseq);
    char *quals = reinterpret_cast<char *>(ref.base + ref.pqual);
    for (int i = 0; i < slen; i++) {
        char base = bases[i];
        char qual = quals[i];
        // get last 3 bits
        char b = base & 0x07;

        const char q20 = '5';
        const char q30 = '?';

        if (qual >= q30) {
            q20bases++;
            q30bases++;
        } else if (qual >= q20) {
            q20bases++;
        }
    }
}

void SeQc::ProducerSeFastqTask(std::string file, rabbit::fq::FastqDataPool *fastqPool,
                               rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    rabbit::fq::FastqFileReader *fqFileReader;
    fqFileReader = new rabbit::fq::FastqFileReader(file, *fastqPool);
    rabbit::int64 n_chunks = 0;
    while (true) {
        rabbit::fq::FastqDataChunk *fqdatachunk;// = new rabbit::fq::FastqDataChunk;
        fqdatachunk = fqFileReader->readNextChunk();
        if (fqdatachunk == NULL) break;
        n_chunks++;
        //std::cout << "readed chunk: " << n_chunks << std::endl;
        dq.Push(n_chunks, fqdatachunk);
    }
    dq.SetCompleted();
    std::cout << "file " << file << " has " << n_chunks << " chunks" << std::endl;
}

void SeQc::PrintRead(neoReference &ref) {
    std::cout << std::string((char *) ref.base + ref.pname, ref.lname) << std::endl;
    std::cout << std::string((char *) ref.base + ref.pseq, ref.lseq) << std::endl;
    std::cout << std::string((char *) ref.base + ref.pstrand, ref.lstrand) << std::endl;
    std::cout << std::string((char *) ref.base + ref.pqual, ref.lqual) << std::endl;
}

void SeQc::ConsumerSeFastqTask(rabbit::fq::FastqDataPool *fastqPool,
                               rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    long line_sum = 0;
    rabbit::int64 id = 0;
//    std::vector<Reference> data;
    rabbit::fq::FastqDataChunk *fqdatachunk;// = new rabbit::fq::FastqDataChunk;
    while (dq.Pop(id, fqdatachunk)) {
        std::vector<neoReference> data;
        //TODO 11000?
//        data.resize(11000);
        int res = rabbit::fq::chunkFormat(fqdatachunk, data, true);
        line_sum += res;
//        for (auto item:data) {
//            StateInfo(item);
//        }
//        std::cout << res << std::endl;
        fastqPool->Release(fqdatachunk);
    }
    std::cout << "line_sum: " << line_sum << std::endl;
//    for (int i = 0; i < 100; i++) {
//        print_read(data[i]);
//        printf("-----------------------------------------------\n");
//    }
    std::cout << "q20bases : " << q20bases << std::endl;
    std::cout << "q20bases : " << q30bases << std::endl;

}

void SeQc::ProcessSeFastq() {
    CmdInfo *cmd_info = SeQc::cmd_info;
    rabbit::fq::FastqDataPool *fastqPool = new rabbit::fq::FastqDataPool(32, 1 << 22);
    //TODO replace this queue
    rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(64, 1);
    //TODO bind ?
    std::thread producer(
            std::bind(&SeQc::ProducerSeFastqTask, this, cmd_info->in_file_name1_, fastqPool, std::ref(queue1)));
    std::thread **threads = new std::thread *[cmd_info->thread_number_];
    for (int t = 0; t < cmd_info->thread_number_; t++) {
        threads[t] = new std::thread(std::bind(&SeQc::ConsumerSeFastqTask, this, fastqPool, std::ref(queue1)));
    }
    producer.join();
    for (int t = 0; t < cmd_info->thread_number_; t++) {
        threads[t]->join();
    }
    delete fastqPool;
    delete[] threads;
}
