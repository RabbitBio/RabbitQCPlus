#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <algorithm>
#include <time.h>
#include <mutex>


inline char complement(char base) {
    switch (base) {
        case 'A':
        case 'a':
            return 'T';
        case 'T':
        case 't':
            return 'A';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        default:
            return 'N';
    }
}

inline bool starts_with(std::string const &value, std::string const &starting) {
    if (starting.size() > value.size()) return false;
    return equal(starting.begin(), starting.end(), value.begin());
}

inline bool ends_with(std::string const &value, std::string const &ending) {
    if (ending.size() > value.size()) return false;
    return equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline std::string trim(const std::string &str) {
    std::string::size_type pos = str.find_first_not_of(' ');
    if (pos == std::string::npos) {
        return std::string("");
    }
    std::string::size_type pos2 = str.find_last_not_of(' ');
    if (pos2 != std::string::npos) {
        return str.substr(pos, pos2 - pos + 1);
    }
    return str.substr(pos);
}

inline int split(const std::string &str, std::vector<std::string> &ret_, std::string sep = ",") {
    if (str.empty()) {
        return 0;
    }

    std::string tmp;
    std::string::size_type pos_begin = str.find_first_not_of(sep);
    std::string::size_type comma_pos = 0;

    while (pos_begin != std::string::npos) {
        comma_pos = str.find(sep, pos_begin);
        if (comma_pos != std::string::npos) {
            tmp = str.substr(pos_begin, comma_pos - pos_begin);
            pos_begin = comma_pos + sep.length();
        } else {
            tmp = str.substr(pos_begin);
            pos_begin = comma_pos;
        }

        ret_.push_back(tmp);
        tmp.clear();
    }
    return 0;
}

inline std::string replace(const std::string &str, const std::string &src, const std::string &dest) {
    std::string ret;

    std::string::size_type pos_begin = 0;
    std::string::size_type pos = str.find(src);
    while (pos != std::string::npos) {
        ret.append(str.data() + pos_begin, pos - pos_begin);
        ret += dest;
        pos_begin = pos + 1;
        pos = str.find(src, pos_begin);
    }
    if (pos_begin < str.length()) {
        ret.append(str.begin() + pos_begin, str.end());
    }
    return ret;
}

inline std::string reverse(const std::string &str) {
    std::string ret(str.length(), 0);
    for (int pos = 0; pos < str.length(); pos++) {
        ret[pos] = str[str.length() - pos - 1];
    }
    return ret;
}

inline std::string basename(const std::string &filename) {
    std::string::size_type pos = filename.find_last_of('/');
    if (pos == std::string::npos)
        return filename;
    else if (pos == filename.length() - 1)
        return "";  // a bad filename
    else
        return filename.substr(pos + 1, filename.length() - pos - 1);
}

inline std::string dirname(const std::string &filename) {
    std::string::size_type pos = filename.find_last_of('/');
    if (pos == std::string::npos) {
        return "./";
    } else
        return filename.substr(0, pos + 1);
}

inline std::string joinpath(const std::string &dirname, const std::string &basename) {
    if (dirname[dirname.length() - 1] == '/') {
        return dirname + basename;
    } else {
        return dirname + "/" + basename;
    }
}

// Check if a string is a file or directory
inline bool file_exists(const std::string &s) {
    bool exists = false;
    if (s.length() > 0) {
        struct stat status;
        int result = stat(s.c_str(), &status);
        if (result == 0) {
            exists = true;
        }
    }
    return exists;
}

// check if a string is a directory
inline bool is_directory(const std::string &path) {
    bool isdir = false;
    struct stat status;
    // visual studion use _S_IFDIR instead of S_IFDIR
    // http://msdn.microsoft.com/en-us/library/14h5k7ff.aspx
#ifdef _MSC_VER
#define S_IFDIR _S_IFDIR
#endif
    stat(path.c_str(), &status);
    if (status.st_mode & S_IFDIR) {
        isdir = true;
    }
    // #endif
    return isdir;
}

inline void check_file_valid(const std::string &s) {
    if (!file_exists(s)) {
        std::cerr << "ERROR: file '" << s << "' doesn't exist, quit now" << std::endl;
        exit(-1);
    }
    if (is_directory(s)) {
        std::cerr << "ERROR: '" << s << "' is a folder, not a file, quit now" << std::endl;
        exit(-1);
    }
}

inline void check_file_writable(const std::string &s) {
    std::string dir = dirname(s);
    if (!file_exists(dir)) {
        std::cerr << "ERROR: '" << dir << " doesn't exist. Create this folder and run this command again." << std::endl;
        exit(-1);
    }
    if (is_directory(s)) {
        std::cerr << "ERROR: '" << s << "' is not a writable file, quit now" << std::endl;
        exit(-1);
    }
}

// Remove non alphabetic characters from a string
inline std::string str_keep_alpha(const std::string &s) {
    std::string new_str;
    for (size_t it = 0; it < s.size(); it++) {
        if (isalpha(s[it])) {
            new_str += s[it];
        }
    }
    return new_str;
}

// Remove invalid sequence characters from a string
inline std::string str_keep_valid_sequence(const std::string &s) {
    std::string new_str;
    for (size_t it = 0; it < s.size(); it++) {
        if (isalpha(s[it]) || s[it] == '-' || s[it] == '*') {
            new_str += s[it];
        }
    }
    return new_str;
}

inline int find_with_right_pos(const std::string &str, const std::string &pattern, int start = 0) {
    int pos = str.find(pattern, start);
    if (pos < 0)
        return -1;
    else
        return pos + pattern.length();
}

inline void str2upper(std::string &s) { transform(s.begin(), s.end(), s.begin(), (int (*)(int))toupper); }

inline void str2lower(std::string &s) { transform(s.begin(), s.end(), s.begin(), (int (*)(int))tolower); }

inline char num2qual(int num) {
    if (num > 127 - 33) num = 127 - 33;
    if (num < 0) num = 0;

    char c = num + 33;
    return c;
}

inline void error_exit(const std::string &msg) {
    std::cerr << "ERROR: " << msg << std::endl;
    exit(-1);
}

extern std::mutex logmtx;
inline void loginfo(const std::string s) {
    logmtx.lock();
    time_t tt = time(NULL);
    tm *t = localtime(&tt);
    std::cerr << "[" << t->tm_hour << ":" << t->tm_min << ":" << t->tm_sec << "] " << s << std::endl;
    logmtx.unlock();
}

#endif /* UTIL_H */
