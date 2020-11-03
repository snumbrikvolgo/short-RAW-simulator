#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <climits>
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
    int64_t given = 0;
    double rate_square;
    double rate;

    int64_t  mean_active;
    int64_t drop;
    int64_t stay;

    double total_delay;
    double square_delay;
    double total_free = 0;

    double free_sum = 0;
    double free_square = 0;
    double delays_sum = 0;
    double delays_sum_square = 0;
    double rates_sum = 0;
    double rates_sum_square = 0;

    int64_t valid_exp = 0;
    int64_t got_exp = 0;


    Counter counter;
};


void PrintAnswer(int N, int K, double lambda, int per_num, int exps_num,
                 SimAnswer const & answer, AccessCategory ac,
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
                   AccessCategory AC, int K,
                   int exps_num, uint64_t seed);


int DropFrame(std::vector<STA> *const stas_ptr, double Tslot, double slot_start_time);


int main() {

    int N = 50;
    int pers_num =  100000, exps_num = 1;
    std::vector <double> ls = {1.00000000e-08};//, 1.26485522e-08, 1.59985872e-08, 2.02358965e-08,
//    2.55954792e-08, 3.23745754e-08, 4.09491506e-08, 5.17947468e-08,
//    6.55128557e-08, 8.28642773e-08, 1.04811313e-07, 1.32571137e-07,
//    1.67683294e-07, 2.12095089e-07, 2.68269580e-07, 3.39322177e-07,
//    4.29193426e-07, 5.42867544e-07, 6.86648845e-07, 8.68511374e-07,
//    1.09854114e-06, 1.38949549e-06, 1.75751062e-06, 2.22299648e-06,
//    2.81176870e-06, 3.55648031e-06, 4.49843267e-06, 5.68986603e-06,
//    7.19685673e-06, 9.10298178e-06, 1.15139540e-05, 1.45634848e-05,
//    1.84206997e-05, 2.32995181e-05, 2.94705170e-05, 3.72759372e-05,
//    4.71486636e-05, 5.96362332e-05, 7.54312006e-05, 9.54095476e-05,
//    1.20679264e-04, 1.52641797e-04, 1.93069773e-04, 2.44205309e-04,
//    3.08884360e-04, 3.90693994e-04, 4.94171336e-04, 6.25055193e-04,
//    7.90604321e-04, 1.00000000e-03
//    };

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
    int K = 15;
    for (int i = 0; i < 1; i++) {
       for (int j = 0; j <= 0; j++) {                                                                                             //ТЕСТЫ ДЛЯ РАЗНЫХ ЛЯМБДА И К
           std::cout << i << " " << j << "\n";
           double Traw = Ts + (K) * Te; //Т.к. в аналит. модели K пустых слотов до успеха
           double Tper = Traw * 10;
           double lambda = ls[i];
           AccessCategory AC(16, 1024, 7, 0); //1 - min, 2 - max, 3 - max RL, AI.. = 0
           auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, K, exps_num, 0U);
//           PrintAnswer(N, K, lambda, pers_num, exps_num, answer.first, AC,false);


           auto pi_end_distr = answer.second.second.getDistribution();
           for ( size_t n = 0; n < pi_end_distr.size(); n++ ) {
               std::cout << pi_end_distr[n] << ",";

           }
            std::cout << answer.first.total_delay / static_cast<float>(answer.first.deliv);
       }
       }
//    for (int i = 0; i < 16; i += 10) {
//        for (int j = 1; j <= 9; j++) {                                                                                             //ТЕСТЫ ДЛЯ РАЗНЫХ ЛЯМБДА И К
//            std::cout << i << " " << j << "\n";
//            double Traw = Ts + (K) * Te; //Т.к. в аналит. модели K пустых слотов до успеха
//            double Tper = Traw * 10;
//            double lambda = ls[i];//-log(1 - q) / 18440;//5.532600e-07;//-log(1 - q) / 18440;//ls[i];//= -log(1 - q) / 18440 ;
//
//
//            AccessCategory AC(16, 1024, j, 0); //1 - min, 2 - max, 3 - max RL, AI.. = 0
//            auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, K, exps_num, 0U);
//            PrintAnswer(N, K, lambda, pers_num, exps_num, answer.first, AC,false);
//        }
//        for (int j = 10; j <= 99; j+=10) {                                                                                             //ТЕСТЫ ДЛЯ РАЗНЫХ ЛЯМБДА И К
//            std::cout << i << " " << j << "\n";
//            double Traw = Ts + (K) * Te; //Т.к. в аналит. модели K пустых слотов до успеха
//            double Tper = Traw * 10;
//            double lambda = ls[i];//-log(1 - q) / 18440;//5.532600e-07;//-log(1 - q) / 18440;//ls[i];//= -log(1 - q) / 18440 ;
//
//
//            AccessCategory AC(16, 1024, j, 0); //1 - min, 2 - max, 3 - max RL, AI.. = 0
//            auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, K, exps_num, 0U);
//            PrintAnswer(N, K, lambda, pers_num, exps_num, answer.first, AC,false);
//        }
//        for (int j = 100; j <= 999; j+=100) {                                                                                             //ТЕСТЫ ДЛЯ РАЗНЫХ ЛЯМБДА И К
//            std::cout << i << " " << j << "\n";
//            double Traw = Ts + (K) * Te; //Т.к. в аналит. модели K пустых слотов до успеха
//            double Tper = Traw * 10;
//            double lambda = ls[i];//-log(1 - q) / 18440;//5.532600e-07;//-log(1 - q) / 18440;//ls[i];//= -log(1 - q) / 18440 ;
//
//
//            AccessCategory AC(16, 1024, j, 0); //1 - min, 2 - max, 3 - max RL, AI.. = 0
//            auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, K, exps_num, 0U);
//            PrintAnswer(N, K, lambda, pers_num, exps_num, answer.first, AC,false);
//        }
//        double Traw = Ts + (K) * Te; //Т.к. в аналит. модели K пустых слотов до успеха
//        double Tper = Traw * 10;
//        AccessCategory AC(16, 1024, INT_MAX, 0); //1 - min, 2 - max, 3 - max RL, AI.. = 0
//        auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, ls[i], AC, K, exps_num, 0U);
//        PrintAnswer(N, K, ls[i], pers_num, exps_num, answer.first, AC,false);
//    }

//    for (int i = 0; i < 9; i++) {
//
//        for (int j = 0; j <= 15; j++)
//        {                                                                                             //ТЕСТЫ ДЛЯ РАЗНЫХ ЛЯМБДА И К
//            std::cout << i << " " << j << "\n";
//            double Traw = Ts + (K) * Te; //Т.к. в аналит. модели K пустых слотов до успеха
//            double Tper = Traw * 10;
//            double lambda = ls[i];//= -log(1 - q) / 18440 ;
//            AccessCategory AC(16, 1024, 7, 0); //1 - min, 2 - max, 3 - max RL, AI.. = 0
//            auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, exps_num, 0U);
//            PrintAnswer(N, K, lambda, pers_num, exps_num, answer.first, AC,false);
//            K++;
//
//        }
//        K = 0;
//        q += 0.1;
//    }
////    std::ofstream active_log;
////    active_log.open("D:\\IITP\\Results\\truth_comparison.txt",  std::ios::out | std::ios::app);
////
////    for (int i = 0; i < max; i++)
////    {
////        std::cout << "cycle" << i << " " << N << "\n";
////        auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, exps_num, 0U);
////        //PrintAnswer(N, K, lambda, pers_num, exps_num, answer.first, false);
////
////        active_log << N << "," << K << "," << lambda << "," << answer.first.mean_active/exps_num << ",";      ПРОВЕРКА СТАЦИОНАРНОГО РАСПРЕДЕЛЕНИЯ

////
////    active_log.close();
    return 0;
}


void PrintAnswer(int N, int K, double lambda, int per, int exps_num,
                 SimAnswer const & answer,AccessCategory ac,
                 bool print_distr=true) {

    std::ofstream active_log;
    active_log.open("D:\\IITP\\Results\\W_diff.txt",  std::ios::out | std::ios::app);
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

            << answer.delays_sum / answer.valid_exp << ","
            << answer.rates_sum / exps_num << ","

            << answer.delays_sum_square / answer.valid_exp - (answer.delays_sum / answer.valid_exp)*(answer.delays_sum / answer.valid_exp) << ","
            << (answer.delays_sum_square / answer.valid_exp - (answer.delays_sum / answer.valid_exp)* (answer.delays_sum / answer.valid_exp)) * answer.valid_exp / (answer.valid_exp - 1) << ","
            << answer.rates_sum_square / exps_num - (answer.rates_sum / exps_num)* (answer.rates_sum / exps_num) << ","
            << (answer.rates_sum_square / exps_num - (answer.rates_sum / exps_num)*(answer.rates_sum / exps_num)) * exps_num / (exps_num - 1) << ","

            << answer.total_free / static_cast<float>(answer.given) << ","
            << answer.free_sum / answer.got_exp << ","
            << answer.free_square / answer.given - (answer.free_sum / answer.given)*(answer.free_sum / answer.given) << ","
            << (answer.free_square / answer.given - (answer.free_sum / answer.given)*(answer.free_sum / answer.given)) * answer.given/ (answer.given - 1) << ","
            << ac.CW_min << ","
            << ac.RL << "\n";

    //answer.delays_sum / exps_num -- выборочное среднее для задержки
    //answer.rates_sum / exps_num -- выборочное среднее для пропускной способности
    //answer.delays_sum_square / exps_num - (answer.delays_sum / exps_num)^2 -- выборочная дисперсия для задержки
    //[answer.delays_sum_square / exps_num - (answer.delays_sum / exps_num)^2] * exps_num / (exps_num - 1) -- несмещённая оценка дисперсии для задержки
    //answer.rates_sum_square / exps_num - (answer.rates_sum / exps_num)^2 -- выборочная дисперсия для пропускной способности
    //[answer.rates_sum_square / exps_num - (answer.rates_sum / exps_num)^2] * exps_num / (exps_num - 1) -- несмещённая оценка дисперсии для пропускной способности


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
                   AccessCategory AC, int K,
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

    std::ofstream active_log;
    active_log.open("D:\\IITP\\Results\\time_dump.txt",  std::ios::out | std::ios::app);


    for (int exp_ind = 0; exp_ind < exps_num; ++exp_ind) { //пошли эксперименты

        std::vector<double> poisson_times;
        poisson_times.reserve(N);       //перед началом эксперимента как такового
        for (int i = 0; i != N; ++i) {
            double poisson_time = -Traw - log(uni_01_dist(mt_gen_eng)) / lambda;
            poisson_times.push_back(poisson_time);      //генерируем вектор рождения пакета во времени по Пуассону
        }

        double total_free =0;
        int64_t given = 0;

        for (int per_ind = 0; per_ind < pers_num; ++per_ind) {
           //N_active = CountActiveStas(&STAs);
            //active_log << N - N_active << ",";

            double current_time = per_ind * Tper + Tper - Traw; //на момент начала слота
            for (int i = 0; i != N; ++i) {
                if ((!STAs[i].ifHasPacket()) && (poisson_times[i] < current_time)) { //если пакета нет и пуассоновское время меньше текущего --> возник пакет в промежутке
                    STAs[i].pushPacket(poisson_times[i]);

                    //std::cout  << "drop " << STAs[i].was_drop << " i " << i << "\n";
                    if (STAs[i].was_drop)
                    {
                        double frtime = poisson_times[i] - STAs[i].lastDropTime();
                        //std::cout  << "frtime " << frtime << "\n";
                        total_free += frtime;
                        answer.total_free += frtime;
                        answer.given ++;
                        given++;
                    }
                    poisson_times[i] = std::numeric_limits<double>::quiet_NaN(); //"удаляем" из очереди присвоения времени
                }
            }

            int N_active = CountActiveStas(&STAs);

            RunRawSlot(&STAs, Traw, current_time, init_freeze);
            //int drop = DropFrame(&STAs, Traw, current_time);                    //ПРОВЕРКА БИНОМИНАЛЬНОГО РАСПРЕДЕЛЕНИЯ
            //active_log << drop << ",";;
            int N_active_after  = CountActiveStas(&STAs);

            if (N_active_after - N_active == -1)
            {
                char st;
                st = 's';
                if (st == 's')
                {
                    std::cout << st ;
                }
            }

//            else if (N_active_after - N_active < -1 || N_active_after > N_active)             ПРОВЕРКА УСПЕХА
//                st = 'b';
//            else st = 'c';

            for (int i = 0; i != N; ++i) {
                if ((!STAs[i].ifHasPacket()) && (std::isnan(poisson_times[i]))) { //нет пакета, но он был, то
                    poisson_times[i] = STAs[i].lastDropTime() - log(uni_01_dist(mt_gen_eng)) / lambda;  //время прихода нового пакета тоже сл.в.
                }
            }

            N_active = CountActiveStas(&STAs); // активные на конец периода
            pi_local.increment(N_active);
            if (per_ind == pers_num - 1)
                answer.mean_active += N_active;
        } /////////////////////////конец всех преиодов

        N_active = CountActiveStas(&STAs); //активные на конец эксперимента
        pi_end.increment(N_active);
        answer.stay += N_active;


        int32_t exp_deliv = 0; // суммарное число кадров переданных за эксперимент
        double exp_sum_delay = 0.0; // суммарная задержка всех кадров за эксперимент

        for (auto &sta : STAs) {
            for (auto const &packet : sta.processed_packets) {
                if (Packet::Status::DELIVERED == packet.status) {
                    ++answer.deliv;
                    answer.total_delay += packet.drop_time - packet.birth_time;
                    answer.square_delay += (packet.drop_time - packet.birth_time) * (packet.drop_time - packet.birth_time) ;

                    ++exp_deliv;
                    exp_sum_delay += packet.drop_time - packet.birth_time;

                } else {
                    ++answer.drop;
                }
                active_log << packet.drop_time << "," << packet.birth_time << "\n";
            }


            answer.counter += sta.counter;
            sta.reset();
        }

        if (given != 0)
        {
            answer.got_exp++; // 1 норм эсперимент по выдаче
            double exp_free = total_free/given; //среднее свободное за эксперимент
            answer.free_sum += exp_free;
            //std::cout << exp_free << "\n";
            answer.free_square += exp_free*exp_free;
        }

        if (exp_deliv != 0)
        {
            answer.valid_exp++;
            double exp_delay = exp_sum_delay / exp_deliv; // средняя задержка за эксперимент
            answer.delays_sum += exp_delay;
            answer.delays_sum_square += exp_delay * exp_delay;
        }

        double exp_rate = exp_deliv / (Tper * pers_num); // пропускная способность за эксперимент

        answer.rates_sum += exp_rate;
        answer.rates_sum_square += exp_rate * exp_rate;
        answer.rate += answer.deliv/Tper/pers_num;
        answer.rate_square += answer.deliv/Tper/pers_num*answer.deliv/Tper/pers_num;
    }//конец 1000 экспериментов

    //answer.delays_sum / exps_num -- выборочное среднее для задержки
    //answer.rates_sum / exps_num -- выборочное среднее для пропускной способности
    //answer.delays_sum_square / exps_num - (answer.delays_sum / exps_num)^2 -- выборочная дисперсия для задержки
    //[answer.delays_sum_square / exps_num - (answer.delays_sum / exps_num)^2] * exps_num / (exps_num - 1) -- несмещённая оценка дисперсии для задержки
    //answer.rates_sum_square / exps_num - (answer.rates_sum / exps_num)^2 -- выборочная дисперсия для пропускной способности
    //[answer.rates_sum_square / exps_num - (answer.rates_sum / exps_num)^2] * exps_num / (exps_num - 1) -- несмещённая оценка дисперсии для пропускной способности

    active_log.close();
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
