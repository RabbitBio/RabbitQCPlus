#ifndef CARE_FILE_HELPERS_HPP
#define CARE_FILE_HELPERS_HPP

#include <cstdio>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <unistd.h>

#include <experimental/filesystem>


//#define FILE_HELPERS_DEBUG

namespace filesys = std::experimental::filesystem;

namespace filehelpers{

__inline__
std::size_t getSizeOfFileBytes(const std::string& filename){
    // https://stackoverflow.com/a/32286531
    
    filesys::path p{filename};
    p = filesys::canonical(p);
    return filesys::file_size(p);
}

__inline__
void renameFileSameMount(const std::string& filename, const std::string& newFilename){
#ifdef FILE_HELPERS_DEBUG        
    std::cerr << "Rename " << filename << " to " << newFilename << "\n";
#endif    
    int res = std::rename(filename.c_str(), newFilename.c_str());
    if(res != 0){
        std::perror("rename");
        assert(res == 0);
    }    
}

__inline__
void copyFile(const std::string& filename, const std::string& newFilename){
#ifdef FILE_HELPERS_DEBUG       
    std::cerr << "Copy " << filename << " to " << newFilename << "\n";
#endif    
    std::ifstream src(filename, std::ios::binary);
    std::ofstream dst(newFilename, std::ios::binary);
    assert(bool(src));
    assert(bool(dst));
    if(src && src.rdbuf()->in_avail() > 0){
        dst << src.rdbuf();
    }
    assert(bool(dst));
}

__inline__
void removeFile(const std::string& filename){
#ifdef FILE_HELPERS_DEBUG   
    std::cerr << "Remove " << filename << "\n";
#endif    
    std::ifstream src(filename);
    assert(bool(src));
    int ret = std::remove(filename.c_str());
    if (ret != 0){
        const std::string errormessage = "Could not remove file " + filename;
        std::perror(errormessage.c_str());
    }  
}

__inline__
std::string makeRandomFile(const std::string& nametemplate){
    std::vector<char> filenamevec(nametemplate.begin(), nametemplate.end());
    filenamevec.push_back('\0');
    int tempfd = mkstemp(filenamevec.data());
    if(tempfd == -1){
        perror("makeRandomFile mkstemp");
        throw std::runtime_error("Cannot create random file with template " + nametemplate);
    }
    close(tempfd);

    return {filenamevec.begin(), filenamevec.end()};
}


__inline__ 
bool fileCanBeOpened(const std::string& filename){
    std::ifstream in(filename);
    return bool(in);
}

__inline__ 
void deleteFiles(const std::vector<std::string>& filenames){
    for (const auto& filename : filenames) {
        removeFile(filename);
    }
}

__inline__
std::uint64_t linecount(const std::string& filename){
	std::uint64_t count = 0;
	std::ifstream is(filename);
	if(is){
		std::string s;
		while(std::getline(is, s)){
			++count;
		}
	}
	return count;
}

__inline__
std::string getFileName(std::string filePath){
    filesys::path path(filePath);
    return path.filename().string();
}

} //namespace filehelpers


#ifdef FILE_HELPERS_DEBUG
#undef FILE_HELPERS_DEBUG
#endif

#endif