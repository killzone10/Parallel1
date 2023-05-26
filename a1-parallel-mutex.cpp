#include <iostream>
#include <thread>
#include <random>
#include <chrono>
#include <complex>
#include <sstream>
#include <string>
#include <fstream>
#include <iomanip>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <future>

#include "a1-helpers.hpp"

template <typename T>
class SafeQ
{
private:
    std::queue<T> q; // no other data structures are allowed
    std::mutex m; // mutex
    std::condition_variable cond; // cond variable

public:
    void push(const T &value)
    {
        std::unique_lock<std::mutex> lock{m}; // locking queue - until pushing value
        q.push(value);
    }
    void push (T && value){
        std::unique_lock<std::mutex> lock{m}; // locking queue - until pushing value
        q.push(std::move(value))

    }

    void pop(T &value)
    {

        std::unique_lock<std::mutex> lock{m}; // locking queue
        
        while (q.empty()){
            cond.wait(lock); // unlocking queue - workers wait there

        }
        value = q.front();
        q.pop();
    }

    void notify(){
        cond.notify_one();

    }
    void notify_all(){ // notifies so the producer can do it from the outside
        cond.notify_all();
    }
    

    std::shared_ptr<T> wait_and_pop(T &value)
    {
        // todo: 
        // in a thread-safe way take the front element
        // and pop it from the queue
        // multiple consumers may be accessing this method
        std::unique_lock<std::mutex> lock{m};

        if ( !q.empty() ){
            std::shared_ptr<T> res_ptr (std::make_shared<T>(q.front()));
            value = q.front();

            q.pop();
            return res_ptr;
        }
        return nullptr;
    }

    size_t size()
    {
        return q.size();
    }

    bool empty()
    {
        return q.empty();
    }
     std::queue<T> get_queue(){
        return q;
     }
};

static constexpr int END_MARKER = std::numeric_limits<int>::min(); // stop criterium

/**
 * To be executed by the master thread
 * The function reads numbers from the file
 * and puts them into the given queue
 * 
 * @param[in] filename
 * @param[inout] q
 * @returns The number of produced items
 * 
*/
int producer(std::string filename, SafeQ<int> &q)
{
    int produced_count = 0;

    // while there are entries in the file
    // put numbers into the queue "q"
    std::ifstream ifs(filename);
    int i = 0;

    while (!ifs.eof()) {
        int num;
        ifs >> num;
        q.push(num);
        i++;
        if (i % 50 == 0){ // notify every nth interaiton
            q.notify();
        }
        produced_count++;
    }
    ifs.close();
    return produced_count;
}

struct WorkerArgs {
    std::mutex mutex; // mutex is there - so the workers dont have race in writing to below vars
    int primes = 0;
    int nonprimes = 0;
    double sum = 0.0;
    int consumed_count = 0;
    std::vector<int> number_counts;
    WorkerArgs() : number_counts(10, 0) {}
};

/**
 * To be executed by worker threads
 * The function removes a number from the queue "q"
 * and does the processing
 * Implement 2 versions with atomic and with mutexes
 * extend as needed
 * 
 * @param[inout] q
 * @param[inout] primes
 * @param[inout] nonprimes
 * @param[inout] mean
 * @param[inout] number_counts
 * 
*/
void worker(SafeQ<int> &q,  WorkerArgs& args, int id)
{
    // implement: use synchronization
    // Note: This part may need some rearranging and rewriting
    // the while loop cannot just check for the size, 
    // it has to now wait until the next element can be popped,
    // or it has to terminate if producer has finished and the queue is empty.

    int consumed_count = 0;
    int primes = 0;
    int nonprimes = 0;
    double sum = 0.0;
    std::vector<int> number_counts(10, 0);

    while (true)
    {
        int num = 0;
        
        q.pop(num);

        if (num == END_MARKER ) {

            // std::cout<<q.size()<<std::endl;
           // q.notify();
            break;
        }

        consumed_count++;
    
        if (kernel(num)) {
            primes++;
        } else {
            nonprimes++;
        }

        number_counts[num % 10]++;
        sum += num;
    }

    std::lock_guard<std::mutex> lock{args.mutex}; // simple lock guard is used there
    args.consumed_count += consumed_count;
    args.primes += primes;
    args.nonprimes += nonprimes;
    args.sum += sum;

    for (int i = 0; i < 10; i++) {
        args.number_counts[i] += number_counts[i];
    }
}

int main(int argc, char **argv)
{
    int num_threads = std::thread::hardware_concurrency();
    std::cout<<num_threads  <<std::endl;
    bool no_exec_times = false, only_exec_times = false; // reporting of time measurements
    std::string filename = "input.txt";
    parse_args(argc, argv, num_threads, filename, no_exec_times, only_exec_times);

    // The actuall code
    // int primes = 0, nonprimes = 0, count = 0;
    // int consumed_count = 0;
    // double mean = 0.0, sum = 0.0;
    // // vector for storing numbers ending with different digits (0-9)
    // std::vector<int> number_counts(10, 0);
    
    // Queue that needs to be made safe 
    // In the simple form it takes integers 
    SafeQ<int> q;
    
    // put you worker threads here
    std::vector<std::thread> workers;

    // time measurement
    auto t1 = std::chrono::high_resolution_clock::now();
    
    // implement: call the producer function with futures/async 
    auto prod = std::async(std::launch::async, producer, filename, std::ref(q));
    WorkerArgs args;

    // implement: spawn worker threads - transform to spawn num_threads threads and store in the "workers" vector
    for (int i = 0; i < num_threads; ++i) {
        workers.emplace_back(worker, std::ref(q), std::ref(args), i);
    }

    int produced_count = prod.get();

    for (int i = 0; i < num_threads; ++i) {
        q.push(END_MARKER);
        q.notify();

    }
    
    for (auto& w : workers) {
        if (w.joinable()) {
            w.join();
        }
    }
    
    double mean = args.sum / args.consumed_count;
    // end time measurement
    auto t2 =  std::chrono::high_resolution_clock::now();

    // do not remove
    if (produced_count != args.consumed_count) {
         std::cout << "[error]: produced_count (" << produced_count << ") != consumed_count (" << args.consumed_count << ")." <<  std::endl;
    }

    // priting the results
    print_output(num_threads, args.primes, args.nonprimes, mean, args.number_counts, t1, t2, only_exec_times, no_exec_times);

    return 0;
}