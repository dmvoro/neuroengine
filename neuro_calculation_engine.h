#ifndef  ___CalculationEngine___
#define  ___CalculationEngine___

//измененная версия класса Олега для нужд обучения нейроклассификатора

#include "ipps.h"
#include "complextypes.h"
#include "frmsdsp.h"

class TFloatIntelFFT;
class TRealFloatIntelFFT;

#include <dsp/global>


struct OldValues  // необходимо для использования в старых алгоритмах вычисления характеристик
{
    unsigned long beg_CSIG;
    unsigned long end_CSIG;
    unsigned long beg_INP;
    int tek_VVSO;  //
    int tek_VVSO1;	// номер последнего считанного блока ДСП25
    int kol_CSIG;
    float fm_sh;      // dima 08.12.

    float sum;  // сумма эспк за 1 цикл
    float apor; // адаптивный порог
    float kpor; // коэффициент apor
    int i, j, err;
    int idxn; /* новый для RAB_PR idx */
    int k;
    int stt;
    int koc; /* номер реал.1-го уровня */
    int kots; /* кол.отсчетов, кол.циклов 2-го уровня */
    int polc, polv, rzn, nspk, nraz, rznk;
    int poln, polk, pols, cntr, cntr1, pln, plk, pft;
    int dras;
    int sdras;
    int pors; /* признак наличия сигнала ors */
    int prh, prh1, crh, crh1, prh4, crh4, crh5;
    int kdr;  /* кол.достоверных реализаций */
    int eso, seso, eso1; /* эмс за 1 реализ.,текущая сумма эмс за 3 цикла */
    int esa;
    int ener;   /* признак непрерыв.энергии */
    int kzap, kzapn; /* коэф.заполнения спк */
    int kzap1, kzat, pnes, knes;
    int kgr, kgrk;
    int nctr, nost1; /* номер строки массива moct */
    int kost1;
    int tax;
    int hax, hdx, hcx, hfd, hfx, hcx1, hcx4, hax1, hdx1;    /* вспомогательные */
    int prob; /* признак обнаружения ВП --------------*/
    int idx;
    int pr_am;
    float correlcoeff;
    int freqmanpr;



    /* МАССИВЫ */
    float  *enrg;   /* сум.enf/ener.выход mmax3.и ps1.---256 скв*/
    float  *enrg1;  /* сум.enf/kots.выход mmax1 и max1.--256 скв*/
    float  *enf;    /* исходный массив СПК(sig) ----------256 скв*/
    float  *eph;    /* массив СПК ofm-свертка**2---------256 скв*/
    float  *eph4;   /* массив СПК ofm4-свертка**4-------256 скв*/
    float  *eno;    /* исх. массив СПК(m1) огиб     ------256 скв*/
    float  *enof;   /* исх. массив СПК(ff1)      ---------256 скв*/
    COMP   *s;      /* массив БПФ ----------------------256 скв*/
    float  *mmax1;  /* исх.enrg1.вых.max1,накап.eso>=0--256 скв*/
    float  *mmax2;  /* вспом.мас нак усреднн СПК   -----256 скв*/
    float  *mmax3;  /* мас.исх-enrg,mmax,mmax2----------256 скв*/
    float  *ff1;    /* исх.мас. мч(sig)           ------ 256 скв*/
    float  *ff2;    /* мас.ff1-УСР2-ПРЖ4       за 12реал.1024 скв*/
    float  *ekf1;   /*                            -------256 ---*/
    float  *ekf2;   /* СПК(ff2) прж.за 12реал.мч.--      512 лк КЧМ и V чм2,4*/
    float  *epf1;   /* СПК(mmax3 на 12-ой рализации) --- 128 лк МКОФ12*/
    float  *epf2;   /* сум.eph(0-255-ОФМ2)(256-511-ОФМ4)-512 скв ОФМ2,4*/
    float  *mvf;    /* мас.ff1-УСР2-УСР4-ПРЖ4 до 300бод-1024 скв*/
    float  *mvf1;   /* вспом.массив для ff1,ff2 -------- 256 лк */
    float  *m1;     /* мас. исходной огибающей ----------256 скв*/
    float  *m2;     /* мас. вспомаг.               ------512 лк ОФМ4 */
    float  *m3;     /* мас. обр. ог(ПРЖ32)         ------256 скв*/
    float  *mmax;    /* мас.сум.enrg/ener СПК------------256 скв*/
    float  *mpss;   /* массив pss -----------------------16 скв*/
    float *HF;	   /* DIM*12-массив частоты  */
    int *peso;      /* массив eso за циклы 1-го уровня----8 скв*/
    int *mson;      /* мас.нак. зон эмс-enf-------------256 скв*/
    int *mson1;     /* мас.зон за 1-реал.---------------256 лк kost,nost*/
    int *mots;      /* мас. отсчетов(эмс) за все реализ. 16 скв*/
    int *moct;      /* 16-массив остатков ------------16 скв*/
    int *moct2;     /* 16-массив стационарный --------16 скв*/
    int *foks;      /* мас.параметров сигнала-----------128 скв*/
    /* 0-  ,1-idx,2-stt,3-RP(rzn),4-CS(polc),5-KP(kgr),6-V(vf2),
    7-eso1(нск\вск),8-pols */
    int *ri;     /* мас.признаков целочислен. --------32 скв*/

    float *rf;  /* мас.признаков с плавующей ---------32 скв*/
    float *mt;	  /* 256*12  накопленной огибающей */
    float *sig;    /* исходный сигнал*/
    float *mv;	 /* мас.ПРЖ4.огиб.за 12реал.--------1024 скв*/
    float *mv1;	 /* вспом. мас. при огибающей  ------256 лк*/

    //============== промежуточные указатели =======================================
    float *uk_f, *uk_f1;
    COMP *uk_c;
    int *uk_i, *uk_i1;
    unsigned int kol_cikl;
    int okno; /* нач.окна спк,кон.окна спк,значение окна */
    float hlp1, hlp2, hlp3, hlp4, max, max1, min, min1, vf2, sh;       //  + sh - D.Z. 06.04
    float temp;	/* для пром.множителей в циклах for() */
    float nost, kost, psps, kgr3, ps1, ps2, psm, pss, pssd, pssa, sump, maxp, minp; /*пер энергия и сумма=pssa+maxp */
    float kgr2, kgrd, srd;
    float poro, pniz, es2;
    int pr_k4m;		/* признак определения КЧМ геродот -----------*/
    int i1, i2, i3;
    float xx, ph;   	 /*   для хранения предыдущей фазы ф-ции CPhasa     */
    float f1;          /*   для хранения предыдущей фазы ф-ции CAmplFreq                */

};

struct SPValues //структура для хранения посчитанных характеристик спектра
{
    int maxOfSp; //максимум спектра
    float maxOfSpVal; //значение максимума спектра
    int secMaxOfSp; //второй максимум спектра
    float secMaxOfSPVal; //значение второго максимума спектра
    int thirdMaxOfSp; //третий максимум спектра
    float thirdMaxOfSpVal; //значение третьего максимума спектра
    int amountOfSamples09; //количество отсчетов спектра выше порога 0.9 от максимума
    int amountOfSamples08; //0.8 от максимума
};



namespace Dsp
{
    class Spectrum;
    class Spectrum_32f;
}


    const int typicalSpectrPower = 10;
    class  NeuroCalculationEngine
    {
    public:
        NeuroCalculationEngine(void);
        ~NeuroCalculationEngine(void);

        //сбрасывает все расчитанные массивы и указатель на исходный сигнал
        void    reset(void);

        //!ПЕРЕД НАЧАЛОМ ЛЮБОЙ ОПЕРАЦИИ (после RESET) НЕЗАБУДЬ СДЕЛАТЬ ОПЕРАЦИЮ setSignal
        //!запоминаем блок сигнала для дальнейшей работы
        void    setSignal(const Ipp32fc *data, int source_size, int sampleRate, int band);
        void    setSignal(NeuroCalculationEngine &ce);
        //!добавить/сдвинуть сигнальный буфер 
        /*! от исходного буфера остается tail_size число отчетов
            добавляются отчеты к этому остатку
            int sampleRate,int band - берутся из ранее вызванной setSignal
            при необходимости выполняется пересчет всех посчитанных ранее величин
                 (в текущей версии относящихся в сигналу envelope,phase,  ) - для экономии быстродействия эти буфера также сдвигаются
            */
        void    append(int tail_size, const Ipp32fc *data, int source_size);
        //!возвращает исходный радиосигнал
        const Ipp32fc   *getSourceSignal(void);
        //!возвращает частоту дискретизации радиосигнала
        int     getSampleRate(void);
        //!возвращает полосу радиосигнала
        int     getBand(void);
        //!возращает размер сигнала
        int     getSourceSignalSize(void);

        //операции над сигналом
        //!сигнал в квадрате
        Ipp32fc*    signalSquare(void);
        //!cсигнал в четвертой степени
        Ipp32fc*    signalFourthDegree(void);
        //!сигнал в восьмой степени
        Ipp32fc*    signalEighthDegree(void);
        //!огибающая радиосигнала
        float*      envelope(void);
        //!фаза радиосигнала
        float*      phase(void);
        //!частота радиосигнала
        float*      frequency(void);
        //!вычисление корреляционной функции от сигнала
        float*      corrSignalAndAmplitude(int corr_size);


        //!логарифм огибающей радиосигнала
        float*      logEnvelope(void);
        //!массив данных доступный для внешнего использования размер равен setSignal(size)*2 = fsource_signal
        float*      tempArray(void);
        //!вычисление корреляционной функции от огибающей // TODO проверить корректность работы корреляционных функций
        float*      corrEnvelope(int corr_data_size);
        float*      modDerEnvelope(void);

        //доп функции
        float*      avFreq(void);
        float*      avModDerFreq(void);

        //!процедура доступа к спектрам
        Ipp32fc*    fftArr(int spectrum_size_order);
        //1амплитудный спектр
        float* ampSpectrFromSignalSource(int spectrum_size_order);
        //!амплитудный спектр сигнала в квадрате
        float* ampSpectrFromSignalSquare(int spectrum_size_order);
        //!амплитудный спектр сигнала в четвертной степени
        float* ampSpectrFromSignalFourthDegree(int spectrum_size_order);
        //!амплитудный спектр сигнала в восьмой степени
        float* ampSpectrFromSignalEighthDegree(int spectrum_size_order);
        //!амплитудный спектр от огибающей сигнала
        float* ampSpectrFromEnvelop(int spectrum_size_order);
        //!амплитудный спектр от автокорреляции огибающей сигнала
        float* ampSpectrFromCorrEnvelope(int spectrum_size_order);
        //!размер спектра
        int    getSpectrumSize(void);
        //!размер спектра от степени 2
        static int     getSpectrumSize(int spectrum_size_order);
        //!размер спектра в степенях двойки
        int    getSpectrumOrder(void);
        //!размер спектра в степенях двойки
        static int    getSpectrumOrder(int spectrum_size);

        //доп функции
        float* ampSpectrFromAmpSpectrFromSignalSource(int spectrum_size_order);
        float* ampSpectrFromAmpSpectrFromSignalSquare(int spectrum_size_order);
        float* ampSpectrFromAmpSpectrFromSignalFourthDegree(int spectrum_size_order);
        float* ampSpectrFromAmpSpectrFromSignalEightDegree(int spectrum_size_order);
        float* ampSpectrFromModDerEnvelope(int spectrum_size_order);  //спектр модуля производной огибающей
        float* ampSpectrFromEnvelopInv(int spectrum_size_order);
        float* powerSpectrum10avFreq(int spectrum_size_order);
        float* powerSpectrum2avModDerFreq(int spectrum_size_order);
        //!простейший частотный детектор //разница между отчетами для которых вычисляется разница 
        /*!задействует temp_array размер результата всегда будет равен getSourceSignalSize()-offset*/
        float* freqDetector(int offset = 1);
        //!более усложненный вариант частотного детектора - пригоден для анализа многочастотных сигналов
        /*!задействует temp_array размер результата всегда будет равен getSourceSignalSize()-1
        */
        float* freqDetector2(void);
        //!третий вариант частотного детектора - предложенный Казаковым Антоном - надо тестить
        //float* freqDetector3(void);


        //ряд сервисных функций
        //!определение ширины сигнала и центра сигнала - возвращается в герцах,  
        /*!max_hole_size - максимальный размер одиночного провала между частями сигнала
               putZero - если true положить сразу 0*/
        int     findSignalWideInHz(float* data, int size, float &max, float &peakSum, float &signalCenter, int max_hole_size = 0, bool putZero = false);
        //!аналогично findSignalWideInHz но ответ в отчетах
        int     findSignalWideInPoints(float* data, int size, float &max, float &peakSum, float &signalCenter, int max_hole_size = 0, bool putZero = false);

        //!нормировка спектра на его общу. сумму
        void    normBySumm(float* data, int size);
        //!медианный фильтр заданного размера
        void    meanWindow(float* data, int size, int windowSize);
        //!медианный фильтр заданного размера
        void    meanWindow(float* source, float *dest, int size, int windowSize);

        //!получение номера спектрального коэффициента от частоты (если частоты за пределами возвращаются наиболее близкие коэффициенты)
        int     getSpectrPos(float frequency);
        //!получение частоты от спектрального коэффициента (-если неверные коэффициенты вернутся наиболее близкие частоты)
        float   getFrequency(float spectrPos);

        //!получение размера хвоста от предыдущего выполнения процедуры append
        int     getTailSignalSize(void);

        //доп функции
        float*      signalKlsFormat(void);
        Ipp32fc*    autoCorr(int spectrum_size_order);
        float*      autoCorrKlsFormat(void);
        void signalToFloat(const Ipp32fc* fSignal, float* resultSignal, int size_of_signal);
        void invertSp(float* nSpectr, int SizeOfSp);
        int findMax(float* Array, int startpoint, int endpoint, float summ, float k); //- поиск максимума с коэффициентом
        int findMaxWithSpace(float* array, int startpoint, int endpoint, int sec_startpoint, int sec_endpoint, float summ, float k, int first_max); // максимума с пропуском
        int findAnotherPeaksWithSpace(float* array, int startpoint, int endpoint, int sec_startpoint, int sec_endpoint, float summ, float k); //поиск  других максимумов в массиве с заданным порогом поиска, возвращается количество максимумов превышающих порог
        bool findSymmetryPeaks(float* array, int startpoint, int endpoint, int sec_startpoint, int sec_endpoint, float summ, float k, int first_max); //ищет 2 симмтричных относительно друг друга пика относительно основного, если есть оба, возвращает true
        float findEnergy(const Ipp32fc* sig, int min_spectrum_size, int startpoint);
        bool  energyTest(const Ipp32fc* sig, int miniSize, float minEnergy, float maxEnergy, int startpoint);
        void findOffsetAndBand(float* spectrum, int &band, int &offset, int spectrum_size_order);
        int findMaxInRegion(float* Array, int startpoint, int endpoint, float summ, float k);
        float* normAmpSpectrum(float* spectrum,int coeff, int spectrum_size_order); //нормирование спектра для того, чтоб подсчет кепстра не отваливался с ошибкой

    protected:
        //последний добавленный сигнал
        Ipp32fc*    fsource_signal;
        int         fsource_signal_size;
        //!хвост сигнала для работы с непрерывными буферами данных
        int         ftail_signal_size;
        int         ftail_max_signal_size;
        int         fsource_max_signal_size;//вводится максимальный размер массива -для уменьшения кол-ва операций удаляния/выделения памяти
        int         fsample_rate;
        int         fband;

        bool        refreshSignalSize(void);
        void        removeArrays(void);

        //для простоты построения системы все сигнальные образцы имеют длительность fsource_signal_size
        Ipp32fc     *fsignal_square, *fsignal_fourth_degree, *fsignal_eighth_degree;
        bool        fsignal_square_computed, fsignal_fourth_degree_computed, fsignal_eighth_degree_computed;
        float*      fenvelope;
        bool        fenvelope_computed;
        float*      fphase;
        bool        fphase_computed;
        float*      ffrequancy;
        float       fprevios_phase;
        bool        ffrequancy_computed;
        float*      flog_envelope;
        bool        flog_envelope_computed;
        float*      fmod_der_envelope;
        bool        fmod_der_envelope_computed;
        float*      ftemp_sig_array;
        Ipp32fc*    ftemp_sig_array_tfc;
        float*      ffreq_detector;
        bool        ffreq_detector_computed;
        Ipp32fc     previos_signal_point;
        float*      fcorr_amp_signal_and_amp;
        bool        fcorr_amp_signal_and_amp_computed;
        
        float*      fsignal_klsformat; // длина 2*fsource_signal_size
        bool        fsignal_klsformat_computed;

        float*      fautocorr_klsformat; // длина 2*fsource_signal_size
        bool        fautocorr_klsformat_computed;
        Ipp32fc*    ffft_arr;
        bool        ffft_arr_computed;


        Dsp::Spectrum* _fft;
        float*      famplitude_spectrum;
        bool        famplitude_spectrum_computed;
        float*      famplitude_spectrum_signal_square;
        bool        famplitude_spectrum_signal_square_computed;
        float*      famplitude_spectrum_signal_fourth_degree;
        bool        famplitude_spectrum_signal_fourth_degree_computed;
        float*      famplitude_spectrum_signal_eighth_degree;
        bool        famplitude_spectrum_signal_eighth_degree_computed;
        float*      famplitude_spectrum_envelope;
        bool        famplitude_spectrum_envelope_computed;
        float*      ftemp_amplitude_spectrum;
        int         fspectrum_size;

        float*      fcorr_data;
        int         fcorr_data_size;
        bool        fcorr_data_computed;
        Ipp8u*      fcorr_internal_buffer;
        int         fcorr_internal_buffer_size;
        float*      famplitude_spectrum_corr_envelope;
        bool        famplitude_spectrum_corr_envelope_computed;

        float*      fautocorr_signal_source;
        bool        fautocorr_signal_source_computed;
        Ipp32fc*    fautocorr_signal_source_tfc;
        bool        fautocorr_signal_source_tfc_computed;

        float*		famplitude_spectrum_envelope_inv;
        bool		famplitude_spectrum_envelope_inv_computed;
        float*      famplitude_spectrum_mod_der_envelope;
        bool        famplitude_spectrum_mod_der_envelope_computed;
        float*      famplitude_spectrum_amplitude_spectrum_signal_source;
        bool        famplitude_spectrum_amplitude_spectrum_signal_source_computed;
        float*      famplitude_spectrum_amplitude_spectrum_signal_fourth_degree;
        bool        famplitude_spectrum_amplitude_spectrum_signal_fourth_degree_computed;
        float*      famplitude_spectrum_amplitude_spectrum_signal_square;
        bool        famplitude_spectrum_amplitude_spectrum_signal_square_computed;
        float*      famplitude_spectrum_amplitude_spectrum_signal_eight_degree;
        bool        famplitude_spectrum_amplitude_spectrum_signal_eight_degree_computed;
       

        float*		fpower_spectrum;
        bool		fpower_spectrum_computed;
        float*		fpower_spectrum_signal_square;
        bool		fpower_spectrum_signal_square_computed;
        float*      fpower_spectrum_signal_fourth_degree;
        bool        fpower_spectrum_signal_fourth_degree_computed;
        float*      fpower_spectrum_signal_eight_degree;
        bool        fpower_spectrum_signal_eight_degree_computed;
        float*		fpower_spectrum_envelope;
        bool		fpower_spectrum_envelope_computed;
        float*      fpower_spectrum_amplitude_spectrum_signal_fourth_degree;
        bool        fpower_spectrum_amplitude_spectrum_signal_fourth_degree_computed;
        float*      fpower_spectrum_amplitude_spectrum_signal_square;
        bool        fpower_spectrum_amplitude_spectrum_signal_square_computed;
        float*      fpower_spectrum_amplitude_spectrum_signal_eight_degree;
        bool        fpower_spectrum_amplitude_spectrum_signal_eight_degree_computed;
        float*      fpower_spectrum_amplitude_spectrum_signal_source;
        bool        fpower_spectrum_amplitude_spectrum_signal_source_computed;

        //еще два дополнительных спектра от частоты для neuro_kls
        //частота - усреднение частоты по 10 отсчетов - спектр мощности
        //спектр мощности модуля производной усредненной по 2 отсчета частоты сигнала
        float*      av_mod_der_freq;
        bool        av_mod_der_freq_computed;
        float*      av_freq;
        bool        av_freq_computed;
        float*      fpower_spectrum_10av_freq;
        bool        fpower_spectrum_10av_freq_computed;
        float*      fpower_spectrum_2av_mod_der_freq;
        bool        fpower_spectrum_2av_mod_der_freq_computed;

        float*      famplitude_norm_spectrum;
        bool        famplitude_norm_spectrum_computed;


        bool        checkSpectrumSize(int size);
        void        removeSpectrumArrays(void);
        float*      compute_amp_spectr(const Ipp32fc* signal, int signal_size, float *spectrum, int pow);
        float*      compute_amp_spectr(const float* signal, int signal_size, float *spectrum, int pow);
        float*		compute_amp_spectr_from_amp_spectr(const float* signal, int signal_size, float *spectrum, int pow);
        //float*      compute_power_spectr(const Ipp32fc* signal, int signal_size, float *spectrum, int pow);
        //float*      compute_power_spectr(const float* signal, int signal_size, float *spectrum, int pow);

        Ipp32fc*    compute_auto_corr_tfc(Ipp32fc* signal, Ipp32fc* autocorr, int nd);
        Ipp32fc*    compute_auto_corr_tfc_sec(Ipp32fc* signal, Ipp32fc* autocorr, int nd);
        Ipp32fc*    compute_fft_arr(Ipp32fc* signal, Ipp32fc* fftarr, int power);

        float*      computeAverage(float* signal, int signal_size, int average_number);

        Dsp::Spectrum_32f* _r_fft;
        TRealFloatIntelFFT* _r_fft_old;
        TFloatIntelFFT* _fft_old;
        OldValues fvalues;
        const int norm_coeff = 1000; //коэффициент нормировки чтоб не падал подсчет кепстра из-за переполнения

    };



#endif