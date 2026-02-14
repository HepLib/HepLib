/**
 * @file 
 * @brief PoolPipe header file
 */
 
#pragma once

#include <vector>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <sys/wait.h>
#include <sys/stat.h>

namespace HepLib {

    class Pool {
    public:
        explicit Pool(const std::vector<void*> & vec);
        void* acquire();
        void release(void* obj);
    private:
        std::vector<void*> ava;
        std::mutex mtx;
        std::condition_variable cv;
    };
    
    class iPool {
    public:
        explicit iPool(int tot);
        explicit iPool();
        int acquire();
        void release(int obj);
        void init(int tot);
    private:
        std::vector<int> ava;
        std::mutex mtx;
        std::condition_variable cv;
        bool inited = false;
    };

    
    class Pipe {
    public:
        explicit Pipe(const std::function<std::string(const std::string &)> & fi);
        ~Pipe();
        std::string run(const std::string &code);
    private:
        int p2c[2], c2p[2];
        pid_t pid;
        std::function<std::string(const std::string &)> f;
    };
    
    class PipePool {
    public:
        explicit PipePool(int size, const std::function<std::string(const std::string &)> & f);
        std::string run(const std::string & code);
        std::vector<std::string> run_all(const std::string &code);
        ~PipePool();
    private:
        Pool* pool;
        std::vector<Pipe*> pipe_ptr_vec;
    };
    
    
}
