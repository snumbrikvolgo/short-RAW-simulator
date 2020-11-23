#ifndef IITP_STATION_H
#define IITP_STATION_H
#pragma once
#include <random>
#include "frame.h"
#include <queue>
#include <cassert>
#include "counter.h"

static int COUNTERS_LIMIT = 5;


struct AccessCategory {
    AccessCategory(int CW_min, int CW_max, int RL, int AIFSN)
            : CW_min(CW_min), CW_max(CW_max), RL(RL), AIFSN(AIFSN) { }

    int CW_min;
    int CW_max;
    int RL;
    int AIFSN;
};


class STA {
public:
    STA(AccessCategory AC, std::mt19937_64 * const gen_ptr);

    void reset();

    void pushPacket(double birth_time);

    bool dropPacket(double drop_time, Packet::Status status);

    double lastDropTime() {
        return last_drop_time;
    }

    void newRaw();

    bool ifTx(int freeze) const {
        return ( !queue.empty() ) && ( freeze >= AC.AIFSN ) && ( 0 == backoff );
    }

    bool ifHasPacket() const {
        return ( !queue.empty() );
    }

    size_t queueLen() const {
        return queue.size();
    }

    int getAIFSN() const {
        return AC.AIFSN;
    }

    void determTxAttempt(bool success, double current_time);

    void edcaSlotBound(int tx_num, int freeze, double current_time);

    std::vector<Packet> processed_packets;
    Counter counter;
    bool was_drop;

    void listenChannel(double W){
        energy_consumption += W;
    }
    double getEnergy() const {
        return  energy_consumption;
    }


private:

    void genBackoff() {
        backoff = std::uniform_int_distribution<int>(0, CW - 1)(*gen_ptr);
        return;
    }

    int retry_counter;
    int CW;
    int backoff;
    std::queue<Packet> queue;

    double last_drop_time;
    double get_time;
    double energy_consumption;


    AccessCategory AC;
    std::mt19937_64 * const gen_ptr;
};

int CountActiveStas(std::vector<STA> * const stas_ptr);

#endif //IITP_STATION_H
