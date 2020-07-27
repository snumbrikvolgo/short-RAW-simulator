#include "station.h"

void STA::reset() {
    processed_packets.clear();
    counter.nullify();
    retry_counter = 0;
    CW = AC.CW_min;
    genBackoff();
    queue = std::queue<Packet>(); //максимальная очередь из пакетов равна 1;
    return;
}

STA::STA(AccessCategory AC, std::mt19937_64 *const gen_ptr)
        : counter({COUNTERS_LIMIT}), retry_counter(0), last_drop_time(0.0), AC(AC), gen_ptr(gen_ptr) {
    CW = AC.CW_min;
    genBackoff();
}

void STA::determTxAttempt(bool success, double current_time) {
    if ( success ) {
        dropPacket(current_time, Packet::Status::DELIVERED);
    } else {
        ++retry_counter;
        if ( AC.RL == retry_counter ) {
            dropPacket(current_time, Packet::Status::DROPPED_RETRY);
        }
    }
    return;
}

void STA::newRaw() {
    counter.increment(queue.size());
    if ( !queue.empty() ) {
        CW = AC.CW_min; //если станция пустая, то какая в целом разница, какой у нее счетчик
        genBackoff();
    }
    return;
}

bool STA::dropPacket(double drop_time, Packet::Status status) {
    if ( queue.empty() ) {
        return false;
    }

    Packet packet = queue.front();
    packet.drop_time = drop_time;
    packet.status = status;
    processed_packets.push_back(packet);

    queue.pop();
    if ( !queue.empty() ) {
        retry_counter = 0;
    }

    last_drop_time = drop_time;

    return true;
}

void STA::pushPacket(double birth_time) {
    if ( queue.empty() ) {
        retry_counter = 0;
    }
    queue.emplace(birth_time);
    return;
}

void STA::edcaSlotBound(int tx_num, int freeze, double current_time) {
    if ( queue.empty() || ( freeze < AC.AIFSN )) {
        return;
    }
    assert(backoff >= 0);

    // this station does not transmit
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
            //CW = std::min(2 * CW, AC.CW_max);
            CW = AC.CW_min; //в этот момент слот завершается, так как была коллизия или успешная передача.
        }
    } else { // success
        dropPacket(current_time, Packet::Status::DELIVERED);
        CW = AC.CW_min;
    }
    genBackoff(); //для дорогих коллизионных и не очень генерируем новые времена отправки
    return;
}


int CountActiveStas(std::vector<STA> * const stas_ptr) {
    int N_active = 0;
    for ( auto const sta : *stas_ptr ) {
        if ( sta.ifHasPacket() ) {
            ++N_active;
        }
    }
    return N_active;
}
