#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include "frame.h"
#include "station.h"
#include "counter.h"

#define Te 52.0 // длина пустого виртуального слота в микросекундах
#define Tc 1064.0 // длина виртуального слота с коллизионной передачей и минимальным AIFS в микросекундах
#define Ts 1064.0 // длина виртуального слота с успешной передачей и минимальным AIFS в микросекундах
#define Tslot_max 246140.0 // максимальный размер слота RAW в микросекундах
#define Data_Length (800.0*1000) // длина полезной части кадра данных 100 байт

struct SimAnswer {
    SimAnswer()
            : deliv(0), drop(0),
              total_delay(0.0),
              counter(COUNTERS_LIMIT) { }

    int64_t deliv;
    int64_t drop;
    double total_delay;
    Counter counter;
};


void PrintAnswer(int N, int exps_num,
                 SimAnswer const & answer,
                 bool print_distr);
int RunRawSlot(std::vector<STA> * const stas_ptr,
               double Tslot, double slot_start_time,
               int init_freeze);
std::pair<SimAnswer, std::pair<Counter, Counter>>
RunSimulation(int N, double Traw, double Tper,
              int pers_num, double lambda,
              AccessCategory AC,
              int exps_num, uint64_t seed);


int main() {
    int N = 50;
    int K = 10;
    double Traw = Ts + (K) * Te; //Т.к. в аналит. модели K пустых слотов до успеха
    double Tper = 100000.0;
    double lambda = 5e-6;
    int pers_num = 1000, exps_num = 1000;
    AccessCategory AC(16, 1024, 7, 0); //1 - min, 2 - max, 3 - max RL, AI.. = 0

    auto answer = RunSimulation(N, Traw, Tper, pers_num, lambda, AC, exps_num, 5489U);

    PrintAnswer(N, exps_num, answer.first, true);

    std::cout << "n\tpi_local\tpi_end\n";
    auto pi_local_distr = answer.second.first.getDistribution();
    auto pi_end_distr = answer.second.second.getDistribution();
    for ( size_t n = 0; n + 1 != pi_local_distr.size(); ++n ) {
        std::cout << n << "\t" << pi_local_distr[n] << "\t" << pi_end_distr[n] << "\n";
    }

    return 0;
}


void PrintAnswer(int N, int exps_num,
                 SimAnswer const & answer,
                 bool print_distr=true) {
    std::cout
            << answer.deliv / static_cast<float>(exps_num * N) << "\t"
            << answer.drop / static_cast<float>(exps_num * N) << "\t"
            << answer.total_delay / static_cast<float>(answer.deliv) << "\n";

    if ( print_distr ) {
        for ( auto const & c : answer.counter.getDistribution() ) {
            std::cout << c << "\t";
        }
        std::cout << "\n";
    }
}

int RunRawSlot(std::vector<STA> *const stas_ptr, double Tslot, double slot_start_time, int init_freeze) {
    for ( auto & sta : *stas_ptr ) {
        sta.newRaw();
    }
    int freeze = init_freeze;
    int e = 0, s = 0, c = 0;

    while ( e * Te + s * Ts + c * Tc <= Tslot - Ts ) {
        int tx_num = 0;
        for ( auto const & sta : *stas_ptr ) {
            if ( sta.ifTx(freeze) ) {
                ++tx_num;
            }
        }

        double current_time = slot_start_time + e * Te + s * Ts + c * Tc;
        if ( 0 == tx_num ) {
            current_time += Te;
        } else if ( 1 == tx_num ) {
            current_time += Ts;
        } else {
            current_time += Tc;
        }

        for ( auto & sta : *stas_ptr ) {
            sta.edcaSlotBound(tx_num, freeze, current_time);
        }

        if ( 0 == tx_num ) {
            ++e;
            ++freeze;
        } else if ( 1 == tx_num ) {
            ++s;
            freeze = 0;
        } else {
            ++c;
            freeze = 0;
        }

        bool smbd_has_packet = false;
        for ( auto const & sta : *stas_ptr ) {
            smbd_has_packet = smbd_has_packet || sta.ifHasPacket();
            if ( smbd_has_packet ) {
                break;
            }
        }
        if ( !smbd_has_packet ) {
            break;
        }
    }
    return s;
}


/**
 * N - число станций, получающих один пакет по событиям Пуассоновского потока
 * Traw - длительность RAW
 * Tper - период RAW
 * pers_num - число периодов в одном эксперименте
 * lambda - интенсивность Пуассоновского потока событий
 * AC - категория доступа, используемая станциями
 * exps_num - число экспериментов
 */
std::pair<SimAnswer, std::pair<Counter, Counter>>
RunSimulation(int N, double Traw, double Tper,
              int pers_num, double lambda,
              AccessCategory AC,
              int exps_num, uint64_t seed=5489U) {

    std::mt19937_64 mt_gen_eng;
    if (0U == seed) {
        mt_gen_eng = std::mt19937_64(clock());
    } else {
        mt_gen_eng = std::mt19937_64(seed);
    }

    std::uniform_real_distribution<double> uni_01_dist;

    std::vector<STA> STAs;
    STAs.reserve(N);
    for (int i = 0; i != N; ++i) {
        STAs.emplace_back(AC, &mt_gen_eng);
    }

    int init_freeze = 0;

    SimAnswer answer;
    Counter pi_local(N);
    Counter pi_end(N);
    int N_active;

    for (int exp_ind = 0; exp_ind < exps_num; ++exp_ind) {
        std::vector<double> poisson_times;
        poisson_times.reserve(N);
        for (int i = 0; i != N; ++i) {
            double poisson_time = -Traw - log(uni_01_dist(mt_gen_eng)) / lambda;
            poisson_times.push_back(poisson_time);
        }

        for (int per_ind = 0; per_ind < pers_num; ++per_ind) {
            double current_time = per_ind * Tper + Tper - Traw;
            for (int i = 0; i != N; ++i) {
                if ((!STAs[i].ifHasPacket()) && (poisson_times[i] < current_time)) {
                    STAs[i].pushPacket(poisson_times[i]);
                    poisson_times[i] = std::numeric_limits<double>::quiet_NaN();
                }
            }
            RunRawSlot(&STAs, Traw, current_time, init_freeze);
            for (int i = 0; i != N; ++i) {
                if ((!STAs[i].ifHasPacket()) && (std::isnan(poisson_times[i]))) {
                    poisson_times[i] = STAs[i].lastDropTime() - log(uni_01_dist(mt_gen_eng)) / lambda;
                }
            }

            N_active = CountActiveStas(&STAs);
            pi_local.increment(N_active);
        }

        N_active = CountActiveStas(&STAs);
        pi_end.increment(N_active);

        for (auto &sta : STAs) {
            for (auto const &packet : sta.processed_packets) {
                if (Packet::Status::DELIVERED == packet.status) {
                    ++answer.deliv;
                    answer.total_delay += packet.drop_time - packet.birth_time;
                } else {
                    ++answer.drop;
                }
            }
            answer.counter += sta.counter;
            sta.reset();
        }
    }

    return std::pair<SimAnswer, std::pair<Counter, Counter>>(answer, std::pair<Counter, Counter>(pi_local, pi_end));
}