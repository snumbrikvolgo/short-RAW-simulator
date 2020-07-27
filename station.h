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

    void edcaSlotBound(int tx_num, int freeze, double current_time) {
        if ( queue.empty() || ( freeze < AC.AIFSN )) {
            return;
        }
        assert(backoff >= 0);

        // this station do not transmit
        if ( backoff > 0 ) {
            --backoff;
            return;
        }

        // this station makes a transmission attempt
        assert(tx_num > 0);
        if ( tx_num > 1 ) { // collision
            ++retry_counter;
            if ( AC.RL == retry_counter ) {
                dropPacket(current_time, Packet::Status::DROPPED_RETRY);
                CW = AC.CW_min;
            } else {
                CW = std::min(2 * CW, AC.CW_max);
            }
        } else { // success
            dropPacket(current_time, Packet::Status::DELIVERED);
            CW = AC.CW_min;
        }
        genBackoff();
        return;
    }

    std::vector<Packet> processed_packets;
    Counter counter;

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

    AccessCategory AC;
    std::mt19937_64 * const gen_ptr;
};

int CountActiveStas(std::vector<STA> * const stas_ptr);

#endif //IITP_STATION_H