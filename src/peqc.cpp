//
// Created by ylf9811 on 2021/7/6.
//
#include <functional>

#include "peqc.h"

PeQc::PeQc(CmdInfo *cmd_info1) {
    cmd_info = cmd_info1;
}

PeQc::~PeQc() {}

void PeQc::StateInfo(neoReference &ref) {

}

void PeQc::PrintRead(neoReference &ref) {

}

void PeQc::ProducerPeFastqTask(std::string file, std::string file2, rabbit::fq::FastqDataPool *fastqPool,
                               rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq) {
    rabbit::fq::FastqFileReader *fqFileReader;
    fqFileReader = new rabbit::fq::FastqFileReader(file, *fastqPool, file2, false);
    int n_chunks = 0;
    int line_sum = 0;
    while (true) {
        rabbit::fq::FastqDataPairChunk *fqdatachunk = new rabbit::fq::FastqDataPairChunk;
        fqdatachunk = fqFileReader->readNextPairChunk();
        if (fqdatachunk == NULL) break;
        n_chunks++;
        //std::cout << "readed chunk: " << n_chunks << std::endl;
        dq.Push(n_chunks, fqdatachunk);
    }

    dq.SetCompleted();
    delete fqFileReader;
    std::cout << "file " << file << " has " << n_chunks << " chunks" << std::endl;
}

void PeQc::ConsumerPeFastqTask(rabbit::fq::FastqDataPool *fastqPool,
                               rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq) {
    long line_sum = 0;
    rabbit::int64 id = 0;
    std::vector<neoReference> data;
    rabbit::fq::FastqDataPairChunk *fqdatachunk = new rabbit::fq::FastqDataPairChunk;
    data.resize(10000);
    while (dq.Pop(id, fqdatachunk)) {
        line_sum += rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk *) (fqdatachunk->left_part), data, true);
        line_sum += rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk *) (fqdatachunk->right_part), data, true);
        fastqPool->Release(fqdatachunk->left_part);
        fastqPool->Release(fqdatachunk->right_part);
    }
}

void PeQc::ProcessPeFastq() {
    CmdInfo *cmd_info = PeQc::cmd_info;
    rabbit::fq::FastqDataPool *fastqPool = new rabbit::fq::FastqDataPool(256, 1 << 22);
    rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> queue1(128, 1);
    std::thread producer(
            std::bind(&PeQc::ProducerPeFastqTask, this, cmd_info->in_file_name1_, cmd_info->in_file_name2_, fastqPool,
                      std::ref(queue1)));
    std::thread **threads = new std::thread *[cmd_info->thread_number_];
    for (int t = 0; t < cmd_info->thread_number_; t++) {
        threads[t] = new std::thread(std::bind(&PeQc::ConsumerPeFastqTask, this, fastqPool, std::ref(queue1)));
    }
    producer.join();

    for (int t = 0; t < cmd_info->thread_number_; t++) {
        threads[t]->join();
    }
    delete fastqPool;
    for (int t = 0; t < cmd_info->thread_number_; t++) {
        delete threads[t];
    }
}
