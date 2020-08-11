#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <tuple>
#include <fstream>
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
            : deliv(0), drop(0),stay(0),
              total_delay(0.0), mean_active(0),
              rate_square(0), rate(0), square_delay(0),
              counter(COUNTERS_LIMIT) { }

    int64_t deliv;
    double rate_square;
    double rate;

    int64_t  mean_active;
    int64_t drop;
    int64_t stay;

    double total_delay;
    double square_delay;
    Counter counter;
};


void PrintAnswer(int N, int K, double lambda, int per_num, int exps_num,
                 SimAnswer const & answer,
                 bool print_distr);
std::tuple<int, int, int> RunRawSlot(std::vector<STA> * const stas_ptr,
               double Tslot, double slot_start_time,
               int init_freeze);

std::pair<SimAnswer, std::pair<Counter, Counter>>
RunSimulation(int N, double Traw, double Tper,
              int pers_num, double lambda,
              AccessCategory AC,
              int exps_num, uint64_t seed);
std::pair<SimAnswer, std::pair<Counter, Counter>>
RunSimulation_Test(int N, double Traw, double Tper,
              int pers_num, double lambda,
              AccessCategory AC,
              int exps_num, uint64_t seed);


int DropFrame(std::vector<STA> *const stas_ptr, double Tslot, double slot_start_time);


int main() {

    int N = 50;
    int max = 50;

//    double Traw = Ts + (K) * Te; //Т.к. в аналит. модели K пустых слотов до успеха
//    double Tper = Traw * 10;
//    double lambda = 7.5e-07;//= -log(1 - q) / 18440 ;
    int pers_num =  100, exps_num = 1000;
    AccessCategory AC(16, 1024, 7, 0); //1 - min, 2 - max, 3 - max RL, AI.. = 0

//    auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, exps_num, 0U);
//
//    PrintAnswer(N, K, lambda,pers_num, exps_num, answer.first, true);
//
//    std::cout << "n\tpi_local\tpi_end\n";
//    auto pi_local_distr = answer.second.first.getDistribution();
//    auto pi_end_distr = answer.second.second.getDistribution();
//    for ( size_t n = 0; n + 1 != pi_local_distr.size(); ++n ) {
//        std::cout << n << "\t" << pi_local_distr[n] << "\t\t" << pi_end_distr[n] << "\n";
//    }




    double q = 0.001;
    int K = 0;
    for (int i = 0; i < 9; i++) {

        for (int j = 0; j <= 15; j++)
        {                                                                                             //ТЕСТЫ ДЛЯ РАЗНЫХ ЛЯМБДА И К
            std::cout << i << " " << j << "\n";
            double Traw = Ts + (K) * Te; //Т.к. в аналит. модели K пустых слотов до успеха
            double Tper = Traw * 10;
            double lambda = -log(1 - q) / 18440;//5.532600e-07;//-log(1 - q) / 18440;//ls[i];//= -log(1 - q) / 18440 ;

            auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, exps_num, 0U);
            PrintAnswer(N, K, lambda, pers_num, exps_num, answer.first, false);

            K++;
        }
        K = 0; //максимум 20 их
        q += 0.001;
    }
    q = 0.01;
    K = 0;
    for (int i = 0; i < 9; i++) {

        for (int j = 0; j <= 15; j++)
        {                                                                                             //ТЕСТЫ ДЛЯ РАЗНЫХ ЛЯМБДА И К
            std::cout << i << " " << j << "\n";
            double Traw = Ts + (K) * Te; //Т.к. в аналит. модели K пустых слотов до успеха
            double Tper = Traw * 10;
            double lambda = -log(1 - q) / 18440;//5.532600e-07;//-log(1 - q) / 18440;//ls[i];//= -log(1 - q) / 18440 ;

            auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, exps_num, 0U);
            PrintAnswer(N, K, lambda, pers_num, exps_num, answer.first, false);

            K++;
        }
        K = 0; //максимум 20 их
        q += 0.01;
    }

    q = 0.1;
    K = 0;
    for (int i = 0; i < 9; i++) {

        for (int j = 0; j <= 15; j++)
        {                                                                                             //ТЕСТЫ ДЛЯ РАЗНЫХ ЛЯМБДА И К
            std::cout << i << " " << j << "\n";
            double Traw = Ts + (K) * Te; //Т.к. в аналит. модели K пустых слотов до успеха
            double Tper = Traw * 10;
            double lambda = -log(1 - q) / 18440;//5.532600e-07;//-log(1 - q) / 18440;//ls[i];//= -log(1 - q) / 18440 ;

            auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, exps_num, 0U);
            PrintAnswer(N, K, lambda, pers_num, exps_num, answer.first, false);

            K++;
        }
        K = 0; //максимум 20 их
        q += 0.1;
    }

        std::vector <double> ls = {5.425710e-08,5.425710e-08
            , 5.425710e-08
            , 7.548410e-08
            , 1.050160e-07
            , 1.085680e-07
            ,  1.461010e-07
            , 1.629340e-07
            ,  2.032600e-07};
    for (int i = 0; i < 9; i++) {

        for (int j = 0; j <= 15; j++)
        {                                                                                             //ТЕСТЫ ДЛЯ РАЗНЫХ ЛЯМБДА И К
            std::cout << i << " " << j << "\n";
            double Traw = Ts + (K) * Te; //Т.к. в аналит. модели K пустых слотов до успеха
            double Tper = Traw * 10;
            double lambda = ls[i];//= -log(1 - q) / 18440 ;

            auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, exps_num, 0U);
            PrintAnswer(N, K, lambda, pers_num, exps_num, answer.first, false);

            K++;
        }
        K = 0; //максимум 20 их
        q += 0.1;
    }



//    std::ofstream active_log;
//    active_log.open("D:\\IITP\\Results\\truth_comparison.txt",  std::ios::out | std::ios::app);
//
//    for (int i = 0; i < max; i++)
//    {
//        std::cout << "cycle" << i << " " << N << "\n";
//        auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, exps_num, 0U);
//        //PrintAnswer(N, K, lambda, pers_num, exps_num, answer.first, false);
//
//        active_log << N << "," << K << "," << lambda << "," << answer.first.mean_active/exps_num << ",";      ПРОВЕРКА СТАЦИОНАРНОГО РАСПРЕДЕЛЕНИЯ
//        auto pi_end_distr = answer.second.second.getDistribution();
//        for ( size_t n = 0; n < pi_end_distr.size(); n++ ) {
//            active_log << pi_end_distr[n] << ",";
//
//        }
//
//        for ( size_t n = pi_end_distr.size(); n <= max; n++ )
//        {
//            active_log << 0 << ",";
//
//        }
//
//        active_log << 0 << "\n";
//
//        N++;
//    }
//
//    active_log.close();
   return 0;
}


void PrintAnswer(int N, int K, double lambda, int per, int exps_num,
                 SimAnswer const & answer,
                 bool print_distr=true) {

    std::ofstream active_log;
    active_log.open("D:\\IITP\\Results\\params_finder_kek.txt",  std::ios::out | std::ios::app);
    active_log
            << answer.deliv / static_cast<float>(exps_num * N) << ","
            << answer.deliv / static_cast<float>(exps_num * per) << ","
            << answer.stay / static_cast<float>(exps_num * N) << ","
            << answer.stay / static_cast<float>(exps_num * per) << ","
            << lambda << ","
            << K << ","
            << answer.drop / static_cast<float>(exps_num * N) << ","
            << answer.drop / static_cast<float>(exps_num * per) << ","
            << answer.total_delay / static_cast<float>(answer.deliv) << ","

            << answer.total_delay / (static_cast<float>(answer.deliv) - 1) << ","
            << answer.square_delay /(static_cast<float>(answer.deliv) - 1) << ","
            << answer.rate  / (static_cast<float>(answer.deliv) - 1)<< ","
            << answer.rate_square /(static_cast<float>(answer.deliv) - 1) << "\n";
//            << answer.deliv << ","
//            << answer.rate/1000<< "\n"; //в секундах


    active_log.close();

    if ( print_distr ) {
        for ( auto const & c : answer.counter.getDistribution() ) {
            std::cout << c << "\t";
        }
        std::cout << "\n";
    }
}

std::tuple<int, int, int> RunRawSlot(std::vector<STA> *const stas_ptr, double Tslot, double slot_start_time, int init_freeze) {
    for ( auto & sta : *stas_ptr ) {
        sta.newRaw();
    }
    int freeze = init_freeze;
    int e = 0, s = 0, c = 0;

    while ( e * Te + s * Ts + c * Tc <= Tslot - Ts ) { //проходим по виртуальным слотам
        int tx_num = 0;
        for ( auto const & sta : *stas_ptr ) {
            if ( sta.ifTx(freeze) ) { //именно ВЫПАДАЯ на этот слот = его начало
                ++tx_num; //сколько передают в этот вирт слот?
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
            sta.edcaSlotBound(tx_num, freeze, current_time); // после этой функции, если была коллизия /успех, слот можно считать завершенным
        }

        if ( 0 == tx_num ) {
            ++e;
            //++freeze;
        } else if ( 1 == tx_num ) {
            ++s;
            //return s;
            //freeze = 0;
        } else {
            ++c;
            //return c;
            //freeze = 0;
        }

        bool smbd_has_packet = false;
        for ( auto const & sta : *stas_ptr ) {
            smbd_has_packet = smbd_has_packet || sta.ifHasPacket();
            if ( smbd_has_packet ) {
                break;              //у кого-то есть пакет
            }
        }
        if ( !smbd_has_packet ) {
            break;                  //если ни у кого его нет , то выходим из функции (для экономии, видимо, времени)
        }
    }
    return std::make_tuple(s,e,c);
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
        STAs.emplace_back(AC, &mt_gen_eng); //создаем N станций со своими счетчиками отсрочки сначала
    }

    int init_freeze = 0;
    std::ofstream active_log;
    active_log.open("D:\\IITP\\Results\\active_log_2.txt",  std::ios::out | std::ios::app );
    SimAnswer answer;
    Counter pi_local(N);        //счетчик активных на конец периода
    Counter pi_end(N);          //счетчик активных на конец эксперимента
    int N_active;

    for (int exp_ind = 0; exp_ind < exps_num; ++exp_ind) { //пошли эксперименты

        std::vector<double> poisson_times;
        poisson_times.reserve(N);       //перед началом эксперимента как такового
        for (int i = 0; i != N; ++i) {
            double poisson_time = -Traw - log(uni_01_dist(mt_gen_eng)) / lambda;
            poisson_times.push_back(poisson_time);      //генерируем вектор рождения пакета во времени по Пуассону
        }

        for (int per_ind = 0; per_ind < pers_num; ++per_ind) {
            double current_time = per_ind * Tper + Tper - Traw; //на момент начала слота
            for (int i = 0; i != N; ++i) {
                if ((!STAs[i].ifHasPacket()) && (poisson_times[i] < current_time)) { //если пакета нет и пуассоновское время меньше текущего --> возник пакет в промежутке
                    STAs[i].pushPacket(poisson_times[i]);
                    poisson_times[i] = std::numeric_limits<double>::quiet_NaN(); //"удаляем" из очереди присвоения времени
                }
            }
//            int given = 0;
//            for (int i = 0; i < N, given!= 6; i++)
//            {
//                if (!STAs[i].ifHasPacket()) {
//                    active_log << "packet given\n";                                                 //РАЗДАЛИ 6 СТАНЦИЯМ
//                    STAs[i].pushPacket(current_time / 2);
//                    given++;
//                }
//            }
//            //МЕНЯЕМ ИМ СЧЕТЧИКИ ОТСРОЧКИ



            N_active = CountActiveStas(&STAs);
            //active_log << N_active << ", ";
            auto per_res = RunRawSlot(&STAs, Traw, current_time, init_freeze);

            int N_active_after = CountActiveStas(&STAs);
//            char st;
//            if (N_active_after - N_active == -1)
//                st = 's';
//            else if (N_active_after - N_active < -1 || N_active_after > N_active)
//                st = 'b';
//            else st = 'c';

            //active_log <<  N_active_after << ", ";


            for (int i = 0; i != N; ++i) {
                if ((!STAs[i].ifHasPacket()) && (std::isnan(poisson_times[i]))) { //нет пакета, но он был, то
                    poisson_times[i] = STAs[i].lastDropTime() - log(uni_01_dist(mt_gen_eng)) / lambda;  //время прихода нового пакета тоже сл.в.
                }
            }

            N_active = CountActiveStas(&STAs); // активные на конец периода
            pi_local.increment(N_active);
        }

        N_active = CountActiveStas(&STAs); //активные на конец эксперимента
        pi_end.increment(N_active);
        for (auto &sta : STAs) {
            for (auto const &packet : sta.processed_packets) {
                if (Packet::Status::DELIVERED == packet.status) {
                    ++answer.deliv;
                    answer.total_delay += packet.drop_time - packet.birth_time;
                    answer.square_delay += (packet.drop_time - packet.birth_time) * (packet.drop_time - packet.birth_time) ;
                } else {
                    ++answer.drop;
                }
            }
            answer.counter += sta.counter;
            //sta.reset();
        }
        answer.rate += answer.deliv/Tper/pers_num;
        answer.rate_square += answer.deliv/Tper/pers_num*answer.deliv/Tper/pers_num;

        for (auto &sta : STAs) {
            sta.reset();
        }
    }
   // active_log.close();

    return std::pair<SimAnswer, std::pair<Counter, Counter>>(answer, std::pair<Counter, Counter>(pi_local, pi_end));
}

std::pair<SimAnswer, std::pair<Counter, Counter>>
RunSimulation_Test(int N, double Traw, double Tper,
              int pers_num, double lambda,
              AccessCategory AC,
              int exps_num, uint64_t seed=0U) {

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
        STAs.emplace_back(AC, &mt_gen_eng); //создаем N станций со своими счетчиками отсрочки сначала
    }

    int init_freeze = 0;

    SimAnswer answer;
    Counter pi_local(N);        //счетчик активных на конец периода
    Counter pi_end(N);          //счетчик активных на конец эксперимента
    int N_active = 0;
    int N_active_after = 0;

//    std::ofstream active_log;
//    active_log.open("D:\\IITP\\Results\\test_null2.txt",  std::ios::out | std::ios::app);


    for (int exp_ind = 0; exp_ind < exps_num; ++exp_ind) { //пошли эксперименты

      std::vector<double> poisson_times;
      poisson_times.reserve(N);       //перед началом эксперимента как такового
      for (int i = 0; i != N; ++i) {
          double poisson_time = -Traw - log(uni_01_dist(mt_gen_eng)) / lambda;
          poisson_times.push_back(poisson_time);      //генерируем вектор рождения пакета во времени по Пуассону
      }

        for (int per_ind = 0; per_ind < pers_num; ++per_ind) {
            N_active = CountActiveStas(&STAs);
            //active_log << N - N_active << ",";

            double current_time = per_ind * Tper + Tper - Traw; //на момент начала слота
            for (int i = 0; i != N; ++i) {
                if ((!STAs[i].ifHasPacket()) && (poisson_times[i] < current_time)) { //если пакета нет и пуассоновское время меньше текущего --> возник пакет в промежутке
                    STAs[i].pushPacket(poisson_times[i]);
                    poisson_times[i] = std::numeric_limits<double>::quiet_NaN(); //"удаляем" из очереди присвоения времени
                }
            }

            N_active = CountActiveStas(&STAs);
            //active_log << N_active << ",";
            RunRawSlot(&STAs, Traw, current_time, init_freeze);
            //int drop = DropFrame(&STAs, Traw, current_time);                    //ПРОВЕРКА БИНОМИНАЛЬНОГО РАСПРЕДЕЛЕНИЯ
            //active_log << drop << ",";;
            N_active_after  = CountActiveStas(&STAs);
            //active_log << N_active_after << "\n";
//            char st;
//            if (N_active_after - N_active == -1)
//                st = 's';
//            else if (N_active_after - N_active < -1 || N_active_after > N_active)             ПРОВЕРКА УСПЕХА
//                st = 'b';
//            else st = 'c';
//
//            active_log <<  N_active_after << ", " << st << "\n";
//

            for (int i = 0; i != N; ++i) {
                if ((!STAs[i].ifHasPacket()) && (std::isnan(poisson_times[i]))) { //нет пакета, но он был, то
                    poisson_times[i] = STAs[i].lastDropTime() - log(uni_01_dist(mt_gen_eng)) / lambda;  //время прихода нового пакета тоже сл.в.
                }
            }

            N_active = CountActiveStas(&STAs); // активные на конец периода
            pi_local.increment(N_active);
            if (per_ind == pers_num - 1)
                answer.mean_active += N_active;
        }

        N_active = CountActiveStas(&STAs); //активные на конец эксперимента
        pi_end.increment(N_active);
        answer.stay += N_active;


        for (auto &sta : STAs) {
            for (auto const &packet : sta.processed_packets) {
                if (Packet::Status::DELIVERED == packet.status) {
                    ++answer.deliv;
                    answer.total_delay += packet.drop_time - packet.birth_time;
                    answer.square_delay += (packet.drop_time - packet.birth_time) * (packet.drop_time - packet.birth_time) ;

                } else {
                    ++answer.drop;
                }
            }


            answer.counter += sta.counter;
            sta.reset();
        }
        //active_log << answer.deliv/Tper/pers_num << "\n";
        answer.rate += answer.deliv/Tper/pers_num;
        answer.rate_square += answer.deliv/Tper/pers_num*answer.deliv/Tper/pers_num;
    }
    //active_log.close();
    return std::pair<SimAnswer, std::pair<Counter, Counter>>(answer, std::pair<Counter, Counter>(pi_local, pi_end));
}

int DropFrame(std::vector<STA> *const stas_ptr, double Tslot, double slot_start_time) {
    for ( auto & sta : *stas_ptr ) {
        sta.newRaw();
    }

    int active = CountActiveStas(stas_ptr);
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    srand((time_t)ts.tv_nsec);
    int dropnum = rand() % (active + 1);
    int dr = dropnum;

    for (int i = 0; i != stas_ptr -> size(), dr > 0; i++)
    {
        if(((*stas_ptr)[i]).ifHasPacket())
        {
            ((*stas_ptr)[i]).dropPacket(slot_start_time, Packet::Status::DELIVERED);
            dr--;
        }
    }

    return dropnum;
}
