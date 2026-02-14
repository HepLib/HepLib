/**
 * @file 
 * @brief PoolPile cpp file
 */

#include "PoolPipe.h"
#include "unistd.h"
#include <stdexcept>
#include <iostream>
#include <arpa/inet.h>

namespace HepLib {

    Pool::Pool(const std::vector<void*> & vec) {
        for (auto item : vec) ava.push_back(item);
    }

    void* Pool::acquire() {
        std::unique_lock<std::mutex> lock(mtx);
        while (ava.empty()) cv.wait(lock);
        void* obj = ava.back();
        ava.pop_back();
        return obj;
    }

    void Pool::release(void* obj) {
        {
            std::lock_guard<std::mutex> lock(mtx);
            ava.push_back(obj);
        }
        cv.notify_one();
    }
    
    iPool::iPool(int tot) {
        for (int i=0; i<tot; i++) ava.push_back(i);
        inited = true;
    }
    
    iPool::iPool() { }
    
    void iPool::init(int tot) {
        for (int i=0; i<tot; i++) ava.push_back(i);
        inited = true;
    }

    int iPool::acquire() {
        if(!inited) throw std::runtime_error("iPool: not initialized.");
        std::unique_lock<std::mutex> lock(mtx);
        while (ava.empty()) cv.wait(lock);
        int obj = ava.back();
        ava.pop_back();
        return obj;
    }

    void iPool::release(int obj) {
        if(!inited) throw std::runtime_error("iPool: not initialized.");
        {
            std::lock_guard<std::mutex> lock(mtx);
            ava.push_back(obj);
        }
        cv.notify_one();
    }

    Pipe::Pipe(const std::function<std::string(const std::string &)> & fi) : f(fi) {
        if (pipe(p2c)==-1 || pipe(c2p)==-1) {
            perror("pipe");
            exit(EXIT_FAILURE);
        }

        pid = fork();
        if (pid == -1) {
            perror("fork");
            exit(EXIT_FAILURE);
        }

        if (pid==0) { // child process
            close(p2c[1]); // close the write side
            close(c2p[0]); // close the read side
            char length_buffer[4];
            uint32_t length;

            while (true) {
                ssize_t bytesRead = read(p2c[0], length_buffer, sizeof(length_buffer));
                if (bytesRead <= 0) break;
                length = ntohl(*(uint32_t*)length_buffer);
                
                char* buffer = new char[length+1];   // +1 for '\0'
                bytesRead = read(p2c[0], buffer, length);
                buffer[bytesRead] = '\0';
                
                std::string cmd(buffer);
                delete[] buffer;
                
                auto res = f(cmd); // cmd -> res

                length = htonl(res.size());
                write(c2p[1], &length, sizeof(length));
                write(c2p[1], res.c_str(), res.size());
            }
            
            close(p2c[0]);
            close(c2p[1]);
            exit(EXIT_SUCCESS);
        } else {
            close(p2c[0]); // close the read
            close(c2p[1]); // close the write
        }
    }

    Pipe::~Pipe() {
        close(p2c[1]); // close the write
        close(c2p[0]); // close the read
        waitpid(pid, NULL, 0);
    }

    std::string Pipe::run(const std::string &code) {
        uint32_t length = htonl(code.size());
        write(p2c[1], &length, sizeof(length));
        write(p2c[1], code.c_str(), code.size());

        char length_buffer[4];
        read(c2p[0], length_buffer, sizeof(length_buffer));
        length = ntohl(*(uint32_t*)length_buffer);

        char* buffer = new char[length+1]; // +1 for '\0'
        ssize_t bytesRead = read(c2p[0], buffer, length);
        buffer[bytesRead] = '\0';
        std::string result(buffer);
        delete[] buffer;

        return result;
    }
    
    PipePool::PipePool(int size, const std::function<std::string(const std::string &)> & f) : pipe_ptr_vec(size) {
        std::vector<void*> void_ptr_vec(size);
        for(int i=0; i<size; i++) {
            pipe_ptr_vec[i] = new Pipe(f);
            void_ptr_vec[i] = pipe_ptr_vec[i];
        }
        pool = new Pool(void_ptr_vec);
    }
    
    PipePool::~PipePool() {
        for(int i=pipe_ptr_vec.size()-1; i>=0; i--) delete pipe_ptr_vec[i];
        delete pool;
    }
    
    std::string PipePool::run(const std::string &code) {
        auto worker = (Pipe*)(pool->acquire());
        auto res = worker->run(code);
        pool->release(worker);
        return res;
    }
    
    std::vector<std::string> PipePool::run_all(const std::string &code) {
        std::vector<std::string> res_vec(pipe_ptr_vec.size());
        for(int i=0; i<res_vec.size(); i++) {
            res_vec[i] = pipe_ptr_vec[i]->run(code);
        }
        return std::move(res_vec);
    }

}
