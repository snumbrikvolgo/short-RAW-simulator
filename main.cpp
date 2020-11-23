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
#define W_idle 2.9 //uJ энергия, потребляемая на прослушивание пустого виртуального слота
#define W_tx 160 // uJ энергия, энергия, потребляемая передающей станицией
#define W_busy 91 //uJ энергия, затрачиваемая на прослушивание занятого виртуального слота

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
    double energy_sum = 0;


    Counter counter;
};

void PrintHeader(std::string const & outfilename);
void PrintAnswer(int N, int K, double lambda, int per_num, int exps_num,
                 SimAnswer const & answer, AccessCategory ac,
                 std::string const & outfilename, bool print_distr);
std::tuple<int, int, int, double> RunRawSlot(std::vector<STA> * const stas_ptr,
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
int main(int argc, char * argv[]) {

    int N = 50;
    int pers_num = 100, exps_num = 100;

    double lambda = std::stod(std::string(argv[1]));
    std::string outfilename = std::string(argv[2]);
    double rl_num = std::stoi(std::string(argv[3]));
    std::cout << rl_num << "\n";
    std::cout << lambda << " " << outfilename << "\n";

    PrintHeader(outfilename);
    for (int j = 2; j <= 64; j++) { // W0
        for (int k = 0; k <= std::min(j - 1, 20); k++) { // K
            std::cout << j << " " << k << "\n";
            double Traw = Ts + (k) * Te; //Т.к. в аналит. модели K пустых слотов до успеха
            double Tper = Traw * 10;
            if (rl_num != 0)
            {
                AccessCategory AC(j, 1024, rl_num, 0); //1 - min, 2 - max, 3 - max RL, AI.. = 0
                auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, k, exps_num, 0U);
                PrintAnswer(N, k, lambda, pers_num, exps_num, answer.first, AC, outfilename, false);
            }

            else
            {
                AccessCategory AC(j, 1024, INT_MAX, 0); //1 - min, 2 - max, 3 - max RL, AI.. = 0
                auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, k, exps_num, 0U);
                PrintAnswer(N, k, lambda, pers_num, exps_num, answer.first, AC, outfilename, false);
            }
//            auto pi_local_distr = answer.second.first.getDistribution();
//            auto pi_end_distr = answer.second.second.getDistribution();
//            for ( size_t n = 0; n + 1 != pi_end_distr.size(); ++n ) {
//                std::cout << n << "\t" << pi_end_distr[n] << "\t\t" << pi_end_distr[n] << "\n";
//            }
        }

    }


    return 0;
}


//int main() {
//
//    int N = 50;
//    int pers_num =  100000, exps_num = 100;
//    std::vector <double> ls = {1.00000000e-08, 1.26485522e-08, 1.59985872e-08, 2.02358965e-08,
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
//
////    auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, exps_num, 0U);
////
////    PrintAnswer(N, K, lambda,pers_num, exps_num, answer.first, true);
////
////    std::cout << "n\tpi_local\tpi_end\n";
////    auto pi_local_distr = answer.second.first.getDistribution();
////    auto pi_end_distr = answer.second.second.getDistribution();
////    for ( size_t n = 0; n + 1 != pi_local_distr.size(); ++n ) {
////        std::cout << n << "\t" << pi_local_distr[n] << "\t\t" << pi_end_distr[n] << "\n";
////    }
//    //int K = 15;
//    //PrintHeader();
//    for (int i = 0; i < 1; i++) {
//        for (int j = 2; j <=
//                         2; j++) {                                                                                             //ТЕСТЫ ДЛЯ РАЗНЫХ ЛЯМБДА И К
//
//            for (int k = 1; k <= 1; k++) {
//                std::cout << i << " " << j << " " << k << "\n";
//                double Traw = Ts + (k) * Te; //Т.к. в аналит. модели K пустых слотов до успеха
//                double Tper = Traw * 10;
//                double lambda = ls[i];
//                AccessCategory AC(j, 1024, 7, 0); //1 - min, 2 - max, 3 - max RL, AI.. = 0
//                auto answer = RunSimulation_Test(N, Traw, Tper, pers_num, lambda, AC, k, exps_num, 0U);
//                //PrintAnswer(N, k, lambda, pers_num, exps_num, answer.first, AC, false);
//            }
//        }
//    }
//
//    return 0;
//}

void PrintHeader(std::string const & outfilename){
    std::ofstream active_log;
    //active_log.open("D:\\IITP\\Results\\energy.txt",  std::ios::out | std::ios::app);
    active_log.open(outfilename,  std::ios::out | std::ios::app);
    active_log << "Delivered per period,Delivered part,Stayed per period,Stayed part,Lambda,K,Dropped per period,Dropped part,Total Delay,Выб ср зад,"<<
               "Выб ср пр,Выб дис зад,Несм дисп зад,Выб дисп пр,Несм дисп пр,Total Free,Выб ср св,Выб дисп св,Несм дисп св,W_0,RL,Power_STA,PLR" << "\n";
    active_log.close();
}
void PrintAnswer(int N, int K, double lambda, int per, int exps_num,
                 SimAnswer const & answer, AccessCategory ac,
                 std::string const & outfilename, bool print_distr) {

    std::ofstream active_log;
    active_log.open(outfilename,  std::ios::out | std::ios::app);
    //active_log.open("D:\\IITP\\Results\\energy.txt",  std::ios::out | std::ios::app);
//    active_log << "Delivered per period,Delivered part,Stayed per period,Stayed part,Lambda,K,Dropped per period,Dropped part,Total Delay,Выб ср зад,"<<
//                  "Выб ср пр,Выб дис зад,Несм дисп зад,Выб дисп пр,Несм дисп пр,Total Free,Выб ср св,Выб дисп св,Несм дисп св,W_0,RL,Power_STA,PLR" << "\n";
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
            << ac.RL << ","
            << answer.energy_sum/static_cast<float>(exps_num * per * N)/static_cast<float>((52*K+1064)*10) << ","
            << answer.drop/(answer.drop + answer.deliv) << "\n";

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

std::tuple<int, int, int, double> RunRawSlot(std::vector<STA> *const stas_ptr, double Tslot, double slot_start_time, int init_freeze) {
    for ( auto & sta : *stas_ptr ) {
        sta.newRaw();
    }
    int freeze = init_freeze;
    int e = 0, s = 0, c = 0;
    double en_raw = 0;

    while ( e * Te + s * Ts + c * Tc <= Tslot - Ts ) { //проходим по виртуальным слотам
        int tx_num = 0;
        for ( auto const & sta : *stas_ptr ) {
            if ( sta.ifTx(freeze) ) { //именно ВЫПАДАЯ на этот слот = его начало
                ++tx_num; //сколько передают в этот вирт слот?
            }
        }
        for ( auto & sta : *stas_ptr ) {
            if (sta.ifHasPacket()){
                if ( !sta.ifTx(freeze) ) { //не передает
                    if (tx_num == 0) // никто вообще не передал
                    {
                        en_raw += W_idle;
                        sta.listenChannel(W_idle); //слушаем пустоту
                    } else {
                        en_raw += W_busy;
                        sta.listenChannel(W_busy); // слушаем занятой канал
                    }
                }
                else {
                    en_raw += W_tx;
                    sta.listenChannel(W_tx);
                }
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
    return std::make_tuple(s,e,c, en_raw);
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
    int N_active_before = 0;
    std::ofstream active_log;
    //active_log.open("D:\\IITP\\Results\\energy.txt",  std::ios::out | std::ios::app);


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

            double current_time = per_ind * Tper + Tper - Traw; //на момент начала слота
            for (int i = 0; i != N; ++i) {
                if ((!STAs[i].ifHasPacket()) && (poisson_times[i] < current_time)) { //если пакета нет и пуассоновское время меньше текущего --> возник пакет в промежутке
                    STAs[i].pushPacket(poisson_times[i]);

                    if (STAs[i].was_drop)
                    {
                        double frtime = poisson_times[i] - STAs[i].lastDropTime();
                        total_free += frtime;
                        answer.total_free += frtime;
                        answer.given ++;
                        given++;
                    }
                    poisson_times[i] = std::numeric_limits<double>::quiet_NaN(); //"удаляем" из очереди присвоения времени
                }
            }
            int N_active_before = CountActiveStas(&STAs);
            auto inf = RunRawSlot(&STAs, Traw, current_time, init_freeze);
            //active_log << N_active << "," << AC.CW_min << "," << K << "," << lambda << "," << std::get<3>(inf) << "\n";

            for (int i = 0; i != N; ++i) {
                if ((!STAs[i].ifHasPacket()) && (std::isnan(poisson_times[i]))) { //нет пакета, но он был, то
                    poisson_times[i] = STAs[i].lastDropTime() - log(uni_01_dist(mt_gen_eng)) / lambda;  //время прихода нового пакета тоже сл.в.
                }
            }

            N_active = CountActiveStas(&STAs); // активные на конец периода

            pi_local.increment(N_active);
            if (per_ind == pers_num - 1)
                answer.mean_active += N_active;

            //answer.energy_sum += std::get<3>(inf);
            //std::cout << per_ind << " " << answer.energy_sum << " "  << std::get<3>(inf) << "\n";
        } /////////////////////////конец всех преиодов
        //std::cout << exp_ind << "\n";
        N_active = CountActiveStas(&STAs); //активные на конец эксперимента
        //active_log << N_active_before << "\n";
        pi_end.increment(N_active_before);
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
            }
            answer.energy_sum += sta.getEnergy();
            answer.counter += sta.counter;
            sta.reset();
        }

        if (given != 0)
        {
            answer.got_exp++; // 1 норм эсперимент по выдаче
            double exp_free = total_free/given; //среднее свободное за эксперимент
            answer.free_sum += exp_free;
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
   // std::cout << answer.energy_sum;

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
