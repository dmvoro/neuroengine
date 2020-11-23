#include "neuro_calculation_engine.h"
#include <math.h>
#include <dsp/spectrum>
#include <fft_library/intel_fft>

using namespace Dsp;

NeuroCalculationEngine::NeuroCalculationEngine(void)
{
    fsource_signal = nullptr;
    fsignal_square = nullptr; 
    fsignal_fourth_degree = nullptr; 
    fsignal_eighth_degree = nullptr;
    fenvelope = nullptr;
    fphase = nullptr;
    ffrequancy = nullptr;
    famplitude_spectrum = nullptr;
    ftemp_amplitude_spectrum = nullptr;
    famplitude_spectrum_signal_square = nullptr;
    famplitude_spectrum_signal_fourth_degree = nullptr;
    famplitude_spectrum_signal_eighth_degree = nullptr;
    famplitude_spectrum_envelope = nullptr;
    famplitude_spectrum_envelope_inv = nullptr;
    famplitude_spectrum_amplitude_spectrum_signal_source = nullptr;
    famplitude_spectrum_amplitude_spectrum_signal_square = nullptr;
    famplitude_spectrum_amplitude_spectrum_signal_fourth_degree = nullptr;
    famplitude_spectrum_amplitude_spectrum_signal_eight_degree = nullptr;
    fpower_spectrum_amplitude_spectrum_signal_source = nullptr;
    fpower_spectrum_amplitude_spectrum_signal_square = nullptr;
    fpower_spectrum_amplitude_spectrum_signal_fourth_degree = nullptr;
    fpower_spectrum_amplitude_spectrum_signal_eight_degree = nullptr;
    famplitude_norm_spectrum = nullptr;
    fspectrum_size = 0;
    flog_envelope = nullptr;
    fmod_der_envelope = nullptr;
    ftemp_sig_array = nullptr;
    fsource_signal_size = 0;
    fsource_max_signal_size = 0;
    fcorr_amp_signal_and_amp = nullptr;

    ftail_signal_size = 0;
    ftail_max_signal_size = 0;

    fcorr_data = nullptr;
    fcorr_data_size = 0;
    fcorr_internal_buffer = nullptr;
    fcorr_internal_buffer_size = 0;
    famplitude_spectrum_corr_envelope = nullptr;
    ffreq_detector = nullptr;

    fsignal_klsformat = nullptr;
    fautocorr_klsformat = nullptr;
    ffft_arr = nullptr;
    fautocorr_signal_source = nullptr;
    fautocorr_signal_source_tfc = nullptr;

    famplitude_spectrum_mod_der_envelope = nullptr;
    famplitude_spectrum_envelope_inv = nullptr;
    fpower_spectrum_10av_freq = nullptr;
    fpower_spectrum_2av_mod_der_freq = nullptr;
    av_freq = nullptr;
    av_mod_der_freq = nullptr;

    fsample_rate = 0;
    fband = 0;

    _fft = new Spectrum;
    _r_fft = new Spectrum_32f;
    reset();
}


NeuroCalculationEngine::~NeuroCalculationEngine(void)
{
    if (_r_fft != nullptr) {
        delete _r_fft;
    }

    if (_fft != nullptr) {
        delete _fft;
    }
    //if (_r_fft_old != nullptr) {
    //    delete _r_fft_old;
    //}
    //if (_fft_old != nullptr) {
    //     delete _fft_old;
    // }

    
    ippsFree(fcorr_internal_buffer);
    removeArrays();
    removeSpectrumArrays();
}

void    NeuroCalculationEngine::removeArrays(void)
{
    if (fsource_signal)
    {
        ippsFree(fsource_signal); fsource_signal = nullptr;
    }
    if (fsignal_square)
    {
        ippsFree(fsignal_square); fsignal_square = nullptr;
    }
    if (fsignal_fourth_degree)
    {
        ippsFree(fsignal_fourth_degree); fsignal_fourth_degree = nullptr;
    }
    if (fsignal_eighth_degree)
    {
        ippsFree(fsignal_eighth_degree); fsignal_eighth_degree = nullptr;
    }
    if (fenvelope)
    {
        ippsFree(fenvelope); fenvelope = nullptr;
    }
    if (fphase)
    {
        ippsFree(fphase); fphase = nullptr;
    }
    if (ffrequancy)
    {
        ippsFree(ffrequancy); ffrequancy = nullptr;
    }
    if (flog_envelope)
    {
        ippsFree(flog_envelope); flog_envelope = nullptr;
    }
    if (fmod_der_envelope)
    {
        ippsFree(fmod_der_envelope); fmod_der_envelope = nullptr;
    }
    if (ftemp_sig_array)
    {
        ippsFree(ftemp_sig_array); ftemp_sig_array = nullptr;
    }
    if (fcorr_data)
    {
        ippsFree(fcorr_data); fcorr_data = nullptr; fcorr_data_size = 0;
    }
    if (ffreq_detector)
    {
        ippsFree(ffreq_detector); ffreq_detector = nullptr;
    }
    if (fcorr_amp_signal_and_amp)
    {
        ippsFree(fcorr_amp_signal_and_amp); fcorr_amp_signal_and_amp = nullptr;
    }
    if (fsignal_klsformat)
    {
        ippsFree(fsignal_klsformat); fsignal_klsformat = nullptr;
    }
    if (fautocorr_klsformat)
    {
        ippsFree(fautocorr_klsformat); fautocorr_klsformat = nullptr;
    }
    if (ffft_arr)
    {
        ippsFree(ffft_arr); ffft_arr = nullptr;
    }
    if (fautocorr_signal_source)
    {
        ippsFree(fautocorr_signal_source); fautocorr_signal_source = nullptr;
    }
    if (fautocorr_signal_source_tfc)
    {
        ippsFree(fautocorr_signal_source_tfc); fautocorr_signal_source_tfc = nullptr;
    }

    
}

void    NeuroCalculationEngine::removeSpectrumArrays(void)
{
    if (famplitude_spectrum)
    {
        ippsFree(famplitude_spectrum); famplitude_spectrum = nullptr;
    }
    if (ftemp_amplitude_spectrum)
    {
        ippsFree(ftemp_amplitude_spectrum); ftemp_amplitude_spectrum = nullptr;
    }

    if (famplitude_spectrum_signal_square)
    {
        ippsFree(famplitude_spectrum_signal_square); famplitude_spectrum_signal_square = nullptr;
    }
    if (famplitude_spectrum_signal_fourth_degree)
    {
        ippsFree(famplitude_spectrum_signal_fourth_degree); famplitude_spectrum_signal_fourth_degree = nullptr;
    }
    if (famplitude_spectrum_signal_eighth_degree)
    {
        ippsFree(famplitude_spectrum_signal_eighth_degree); famplitude_spectrum_signal_eighth_degree = nullptr;
    }
    if (famplitude_spectrum_envelope)
    {
        ippsFree(famplitude_spectrum_envelope); famplitude_spectrum_envelope = nullptr;
    }
    if (famplitude_spectrum_corr_envelope)
    {
        ippsFree(famplitude_spectrum_corr_envelope); famplitude_spectrum_corr_envelope = nullptr;
    }
    
    if (famplitude_spectrum_mod_der_envelope)
    {
        ippsFree(famplitude_spectrum_mod_der_envelope); famplitude_spectrum_mod_der_envelope = nullptr;
    }
    if (famplitude_spectrum_envelope_inv)
    {
        ippsFree(famplitude_spectrum_envelope_inv);  famplitude_spectrum_envelope_inv = nullptr;
    }
    if (fpower_spectrum_10av_freq)
    {
        ippsFree(fpower_spectrum_10av_freq); fpower_spectrum_10av_freq = nullptr;
    }
    if (fpower_spectrum_2av_mod_der_freq)
    {
        ippsFree(fpower_spectrum_2av_mod_der_freq); fpower_spectrum_2av_mod_der_freq = nullptr;
    }
    if (famplitude_spectrum_amplitude_spectrum_signal_square)
    {
        ippsFree(famplitude_spectrum_amplitude_spectrum_signal_square); famplitude_spectrum_amplitude_spectrum_signal_square = nullptr;
    }
    if (famplitude_spectrum_amplitude_spectrum_signal_fourth_degree)
    {
        ippsFree(famplitude_spectrum_amplitude_spectrum_signal_fourth_degree); famplitude_spectrum_amplitude_spectrum_signal_fourth_degree = nullptr;
    }
    if (famplitude_spectrum_amplitude_spectrum_signal_eight_degree)
    {
        ippsFree(famplitude_spectrum_amplitude_spectrum_signal_eight_degree); famplitude_spectrum_amplitude_spectrum_signal_eight_degree = nullptr;
    }
    if (famplitude_spectrum_amplitude_spectrum_signal_source)
    {
        ippsFree(famplitude_spectrum_amplitude_spectrum_signal_source); famplitude_spectrum_amplitude_spectrum_signal_source = nullptr;
    }
    if (famplitude_norm_spectrum)
    {
        ippsFree(famplitude_norm_spectrum); famplitude_norm_spectrum = nullptr;
    }
}

void    NeuroCalculationEngine::reset(void)
{   //сброс вычисленных массивов сигналов
    fsignal_square_computed = false; fsignal_fourth_degree_computed = false; fsignal_eighth_degree_computed = false;
    fenvelope_computed = false;
    fphase_computed = false;
    flog_envelope_computed = false;
    fcorr_data_computed = false;
    ffrequancy_computed = false;
    fprevios_phase = 0;
    ffreq_detector_computed = false;
    fcorr_amp_signal_and_amp_computed = false;
    fsignal_klsformat_computed = false;
    fautocorr_klsformat_computed = false;
    ffft_arr_computed = false;
    fautocorr_signal_source_computed = false;
    fautocorr_signal_source_tfc_computed = false;
    fmod_der_envelope_computed = false;

    //сброс вычисленных массивов спектров
    famplitude_spectrum_computed = false;
    famplitude_spectrum_signal_square_computed = false;
    famplitude_spectrum_envelope_computed = false;
    famplitude_spectrum_signal_fourth_degree_computed = false;
    famplitude_spectrum_signal_eighth_degree_computed = false;
    famplitude_spectrum_corr_envelope_computed = false;

    famplitude_spectrum_amplitude_spectrum_signal_source_computed = false;
    famplitude_spectrum_amplitude_spectrum_signal_square_computed = false;
    famplitude_spectrum_amplitude_spectrum_signal_fourth_degree_computed = false;
    famplitude_spectrum_amplitude_spectrum_signal_eight_degree_computed = false;
    fpower_spectrum_amplitude_spectrum_signal_source_computed = false;
    fpower_spectrum_amplitude_spectrum_signal_square_computed = false;
    fpower_spectrum_amplitude_spectrum_signal_fourth_degree_computed = false;
    fpower_spectrum_amplitude_spectrum_signal_eight_degree = false;

    famplitude_spectrum_mod_der_envelope_computed = false;
    famplitude_spectrum_envelope_inv_computed = false;
    fpower_spectrum_10av_freq_computed = false;
    fpower_spectrum_2av_mod_der_freq_computed = false;

    av_freq_computed = false;
    av_mod_der_freq_computed = false;
    famplitude_norm_spectrum_computed = false;

    //сброс данных о рассчитанном хвосте
    ftail_signal_size = 0;
    //забываем о сигнале
    fsource_signal_size = 0;
    //сбрасываем предыдущий отчет
    previos_signal_point.re = 0.0f;
    previos_signal_point.im = 0.0f;
}

//запоминаем блок сигнала для дальнейшей работы
void    NeuroCalculationEngine::setSignal(const Ipp32fc *data, int size, int sampleRate, int band)
{
    reset();
    if (fsource_max_signal_size < size)
    {
        fsource_max_signal_size = size;
        refreshSignalSize();
        fsource_signal = ippsMalloc_32fc(fsource_max_signal_size);
    }
    fsource_signal_size = size;
    fsample_rate = sampleRate;
    fband = band;
    ippsCopy_32fc(data, fsource_signal, size);
}

void    NeuroCalculationEngine::setSignal(NeuroCalculationEngine &ce)
{
    setSignal(ce.getSourceSignal(), ce.getSourceSignalSize(), ce.getSampleRate(), ce.getBand());
}


void    NeuroCalculationEngine::append(int tail_size, const Ipp32fc *data, int source_size)
{
    //копируем хвост
    int size = fsource_signal_size + ftail_signal_size;
    int copy_size = size > tail_size ? tail_size : size;
    int copy_pointer = size - copy_size;
    ftail_signal_size = tail_size;

    //сохраняем указатели на все сигнальные массивы данных для использования в дальнейшем
    Ipp32fc     *fold_source_signal = fsource_signal;
    Ipp32fc     *fold_signal_square = fsignal_square;
    Ipp32fc     *fold_signal_fourth_degree = fsignal_fourth_degree;
    Ipp32fc     *fold_signal_eighth_degree = fsignal_eighth_degree;
    float       *fold_envelope = fenvelope;
    //float       *fold_fmod_der_envelope = fmod_der_envelope
    float       *fold_phase = fphase;
    float       *fold_frequancy = ffrequancy;
    float       *fold_log_envelope = flog_envelope;
    float       *fold_freq_detector = ffreq_detector;
    //необходимо отмасштабировать сигнальные массивы данных
    if (fsource_max_signal_size < source_size + ftail_signal_size)
    {
        fsource_max_signal_size = source_size + ftail_signal_size;
        fsource_signal = ippsMalloc_32fc(fsource_max_signal_size);
        //создание новых массивов данных
        if (fsignal_square_computed)
            fsignal_square = ippsMalloc_32fc(fsource_max_signal_size);
        if (fsignal_fourth_degree_computed)
            fsignal_fourth_degree = ippsMalloc_32fc(fsource_max_signal_size);
        if (fsignal_eighth_degree_computed)
            fsignal_eighth_degree = ippsMalloc_32fc(fsource_max_signal_size);
        if (fenvelope_computed)
            fenvelope = ippsMalloc_32f(fsource_max_signal_size);
        //if(fmod_der_envelope_computed)
        //    fmod_der_envelope = ippsMalloc_32f(fsource_max_signal_size);
        if (fphase_computed)
            fphase = ippsMalloc_32f(fsource_max_signal_size);
        if (ffrequancy_computed)
            ffrequancy = ippsMalloc_32f(fsource_max_signal_size);
        if (flog_envelope_computed)
            flog_envelope = ippsMalloc_32f(fsource_max_signal_size);
        if (ffreq_detector_computed)
            ffreq_detector = ippsMalloc_32f(fsource_max_signal_size);
        if (ftemp_sig_array)
        {
            ippsFree(ftemp_sig_array); ftemp_sig_array = nullptr;
        }
    }
    //копируем уже посчитанные блоки данных и удаляем старый буфер по мере необходимости
    ippsMove_32fc(&fold_source_signal[copy_pointer], fsource_signal, copy_size);
    if (fold_source_signal != fsource_signal)
        ippsFree(fold_source_signal);

    if (fsignal_square_computed)
    {
        ippsMove_32fc(&fold_signal_square[copy_pointer], fsignal_square, copy_size);
        if (fold_signal_square != fsignal_square)
            ippsFree(fold_signal_square);
    }
    if (fsignal_fourth_degree_computed)
    {
        ippsMove_32fc(&fold_signal_fourth_degree[copy_pointer], fsignal_fourth_degree, copy_size);
        if (fold_signal_fourth_degree != fsignal_fourth_degree)
            ippsFree(fold_signal_fourth_degree);
    }
    if (fsignal_eighth_degree_computed)
    {
        ippsMove_32fc(&fold_signal_eighth_degree[copy_pointer], fsignal_eighth_degree, copy_size);
        if (fold_signal_eighth_degree != fsignal_eighth_degree)
            ippsFree(fold_signal_eighth_degree);
    }
    if (fenvelope_computed)
    {
        ippsMove_32f(&fold_envelope[copy_pointer], fenvelope, copy_size);
        if (fold_envelope != fenvelope)
            ippsFree(fold_envelope);
    }
    /*if (fmod_der_envelope_computed)
    {
        ippsMove_32f(&fold_fmod_der_envelope[copy_pointer], fmod_der_envelope, copy_size);
        if (fold_envelope != fenvelope)
            ippsFree(fold_envelope);
    }*/
    if (fphase_computed)
    {
        ippsMove_32f(&fold_phase[copy_pointer], fphase, copy_size);
        if (fold_phase != fphase)
            ippsFree(fold_phase);
    }
    if (ffrequancy_computed)
    {
        ippsMove_32f(&fold_frequancy[copy_pointer], ffrequancy, copy_size);
        if (fold_frequancy != ffrequancy)
            ippsFree(fold_frequancy);
    }
    if (flog_envelope_computed)
    {
        ippsMove_32f(&fold_log_envelope[copy_pointer], flog_envelope, copy_size);
        if (fold_log_envelope != flog_envelope)
            ippsFree(fold_log_envelope);
    }
    if (ffreq_detector_computed)
    {
        ippsMove_32f(&fold_freq_detector[copy_pointer], ffreq_detector, copy_size);
        if (fold_freq_detector != ffreq_detector)
            ippsFree(fold_freq_detector);
    }

    //копирование нового блока сигнала
    ippsCopy_32fc(data, &fsource_signal[copy_size], source_size);//fsource_signal_size);

    //int tmp = ftail_signal_size;
    reset();
    //операции непосредственно с сигналом
    fsource_signal_size = source_size;
    ftail_signal_size = copy_size;

}

float*      NeuroCalculationEngine::tempArray(void)
{
    if (!ftemp_sig_array)
        ftemp_sig_array = ippsMalloc_32f(fsource_max_signal_size * 2);
    return ftemp_sig_array;
};


bool    NeuroCalculationEngine::refreshSignalSize(void)
{
    if (!fsource_max_signal_size)
        return false;
    removeArrays();
    //массивы выделяются по мере необходимости  
    return true;
}

Ipp32fc* NeuroCalculationEngine::signalSquare(void)
{
    if (!fsource_signal)
        return nullptr;
    if (!fsignal_square_computed)
    {
        if (!fsignal_square)
            fsignal_square = ippsMalloc_32fc(fsource_max_signal_size);
        ippsMul_32fc(fsource_signal, fsource_signal, &fsignal_square[ftail_signal_size], fsource_signal_size);
        fsignal_square_computed = true;
    }
    return fsignal_square;
}

Ipp32fc* NeuroCalculationEngine::signalFourthDegree(void)
{
    if (!fsource_signal)
        return nullptr;
    if (!fsignal_fourth_degree_computed)
    {
        if (!fsignal_fourth_degree)
            fsignal_fourth_degree = ippsMalloc_32fc(fsource_max_signal_size);
        ippsMul_32fc(signalSquare(), signalSquare(), &fsignal_fourth_degree[ftail_signal_size], fsource_signal_size);
        fsignal_fourth_degree_computed = true;
    }
    return fsignal_fourth_degree;
}

Ipp32fc* NeuroCalculationEngine::signalEighthDegree(void)
{
    if (!fsource_signal)
        return nullptr;
    if (!fsignal_eighth_degree_computed)
    {
        if (!fsignal_eighth_degree)
            fsignal_eighth_degree = ippsMalloc_32fc(fsource_max_signal_size);
        ippsMul_32fc(signalFourthDegree(), signalFourthDegree(), &fsignal_eighth_degree[ftail_signal_size], fsource_signal_size);
        fsignal_eighth_degree_computed = true;
    }
    return fsignal_eighth_degree;
}

float*  NeuroCalculationEngine::envelope(void)
{
    if (!fsource_signal)
        return nullptr;
    if (!fenvelope_computed)
    {
        if (!fenvelope)
            fenvelope = ippsMalloc_32f(fsource_max_signal_size);
        ippsMagnitude_32fc(fsource_signal, &fenvelope[ftail_signal_size], fsource_signal_size);
        fenvelope_computed = true;
    }
    return fenvelope;
}

float*  NeuroCalculationEngine::phase(void)
{
    if (!fsource_signal)
        return nullptr;
    if (!fphase_computed)
    {
        if (!fphase)
            fphase = ippsMalloc_32f(fsource_max_signal_size);
        ippsPhase_32fc(fsource_signal, &fphase[ftail_signal_size], fsource_signal_size);
        fphase_computed = true;
    }
    return fphase;
}
float*  NeuroCalculationEngine::frequency(void)
{
    if (!fsource_signal)
    {
        return nullptr;
    }
    if (!ffrequancy_computed)
    {
        if (!ffrequancy)
            ffrequancy = ippsMalloc_32f(fsource_max_signal_size);
        float *ph = phase();
        for (int i = 0; i < fsource_signal_size; i++)//можно считать меньше OLBER BUG
        {
            float Freq = (ph[i] - fprevios_phase);
            constexpr float F_IPP_PI = static_cast<float>(IPP_PI);//
            if (Freq > F_IPP_PI)
                Freq -= F_IPP_PI * 2.0f;
            if (Freq <= -F_IPP_PI)
                Freq += F_IPP_PI * 2.0f;
            fprevios_phase = ph[i];
            ffrequancy[i + ftail_signal_size] = Freq;
        }
        ffrequancy_computed = true;
    }
    return ffrequancy;
}

//!вычисление корреляционной функции от сигнала
float*  NeuroCalculationEngine::corrSignalAndAmplitude(int corr_size)
{
    if (corr_size >= getSourceSignalSize())
        return nullptr;
    if (!fcorr_amp_signal_and_amp_computed)
    {
        if (!fcorr_amp_signal_and_amp)
            fcorr_amp_signal_and_amp = ippsMalloc_32f(fsource_max_signal_size);
        if (fcorr_internal_buffer_size < corr_size*static_cast<int>(sizeof(Ipp32fc)))
        {
            ippsFree(fcorr_internal_buffer);
            fcorr_internal_buffer = ippsMalloc_8u_L(corr_size*static_cast<int>(sizeof(Ipp32fc)));
        }
        int steps_count = getSourceSignalSize() - corr_size;
        const Ipp32fc* sig = getSourceSignal();
        Ipp32fc  *fcorr_sig = reinterpret_cast<Ipp32fc*>(fcorr_internal_buffer);
        Ipp32fc  *ftemp_sig = reinterpret_cast<Ipp32fc*>(tempArray());
        //делаем комплексно-сопряженное - TODO надо сделать лучший выбор образца
        float max; int index;
        ippsMaxIndx_32f(reinterpret_cast<const float*>(sig), getSourceSignalSize() * 2, &max, &index);
        index = index - corr_size / 2;
        index = int(index / corr_size)*corr_size;
        if (index > getSourceSignalSize() - corr_size)
            index = getSourceSignalSize() - corr_size;
        for (int i = 0; i < corr_size; i++)
        {
            fcorr_sig[i].re = sig[i + index].re;
            fcorr_sig[i].im = -sig[i + index].im;
        }
        //предполагаем что все слоты активны
        Ipp32fc sum;
        for (int i = 0; i < steps_count; i++)
        {
            ippsMul_32fc(&sig[i], fcorr_sig, ftemp_sig, corr_size);
            ippsSum_32fc(ftemp_sig, corr_size, &sum, ippAlgHintFast);
            fcorr_amp_signal_and_amp[i] = sum.re*sum.re + sum.im*sum.im;
        }
    }
    return fcorr_amp_signal_and_amp;
}

float*  NeuroCalculationEngine::logEnvelope(void)
{
    if (!fsource_signal)
        return nullptr;
    if (!flog_envelope_computed)
    {
        if (!flog_envelope)
            flog_envelope = ippsMalloc_32f(fsource_max_signal_size);
        ippsAddC_32f_I(0.0000000001f, &envelope()[ftail_signal_size], fsource_signal_size);
        ippsLn_32f(&envelope()[ftail_signal_size], &flog_envelope[ftail_signal_size], fsource_signal_size);
        flog_envelope_computed = true;
    }
    return flog_envelope;
}

//вычисление корреляционной функции от огибающей
float*      NeuroCalculationEngine::corrEnvelope(int corr_data_size)
{

    if (!fcorr_data_computed)
    {
        int signal_size = getSourceSignalSize();
        if (corr_data_size != fcorr_data_size)
        {
            ippsFree(fcorr_data);
            fcorr_data_size = corr_data_size;
            fcorr_data = ippsMalloc_32f(fcorr_data_size);

        }
        int size;
        IppStatus st = ippsAutoCorrNormGetBufferSize(signal_size, fcorr_data_size, ipp32f, ippiNormNone, &size);
        Q_ASSERT(ippStsNoErr == st);
        if (fcorr_internal_buffer_size < size)
        {
            ippsFree(fcorr_internal_buffer);
            fcorr_internal_buffer = ippsMalloc_8u_L(size);
        }
        st = ippsAutoCorrNorm_32f(envelope(), signal_size, fcorr_data, fcorr_data_size, ippiNormNone, fcorr_internal_buffer);
        Q_ASSERT(ippStsNoErr == st);

        //выделим производную
        for (int i = 0; i < fcorr_data_size - 1; i++)
        {
            fcorr_data[i] -= fcorr_data[i + 1];
        }
        fcorr_data[fcorr_data_size - 1] = fcorr_data[fcorr_data_size - 2];
        fcorr_data[0] = fcorr_data[1];

        /*
        float *env = envelope();
        for (int i=0;i<corr_data_size;i++)
        {   float sum   = nullptr;
            for (int j=0;j<signal_size-corr_data_size;j++)
            {   sum+=env[j]*env[i+j];
            }
            fcorr_data[i] = sum;
        }
        */
        fcorr_data_computed = true;
    }


    return fcorr_data;
}

bool   NeuroCalculationEngine::checkSpectrumSize(int size)
{
    if (size != fspectrum_size)
    {
        if (!size)
            return false;
        fspectrum_size = size;
        removeSpectrumArrays();
        //TO DO optimization need !!!!
    }
    return true;
}

float*   NeuroCalculationEngine::compute_amp_spectr(const Ipp32fc* signal, int signal_size, float *spectrum, int power)
{
    int spectrum_size = static_cast<int>(pow(2.0f, power));
    if (fspectrum_size > signal_size)
        return nullptr;
    if (!ftemp_amplitude_spectrum)
        ftemp_amplitude_spectrum = ippsMalloc_32f(fspectrum_size);
    _fft->magnutude(signal, spectrum_size, spectrum);
    for (int i = fspectrum_size; i <= signal_size - fspectrum_size; i += fspectrum_size / 2)
    {
        _fft->magnutude(&signal[i], spectrum_size, ftemp_amplitude_spectrum);
        //normBySumm(ftemp_amplitude_spectrum,fspectrum_size);
        ippsAdd_32f_I(ftemp_amplitude_spectrum, spectrum, fspectrum_size);
    }
    normBySumm(spectrum, fspectrum_size);
    return spectrum;
}

float*  NeuroCalculationEngine::compute_amp_spectr(const float* signal, int signal_size, float *spectrum, int power)
{
    int spectrum_size = static_cast<int>(pow(2.0f, power));
    if (fspectrum_size > signal_size)
        return nullptr;
    if (!ftemp_amplitude_spectrum)
        ftemp_amplitude_spectrum = ippsMalloc_32f(fspectrum_size);
    _r_fft->magnutude(signal, spectrum_size, spectrum);
    for (int i = fspectrum_size; i < signal_size - fspectrum_size; i += fspectrum_size)
    {
        _r_fft->magnutude(&signal[i], spectrum_size, ftemp_amplitude_spectrum);
        ippsAdd_32f_I(ftemp_amplitude_spectrum, spectrum, fspectrum_size);
    }
    normBySumm(spectrum, fspectrum_size);
    return spectrum;
}
//#include "qfile.h"
float* NeuroCalculationEngine::ampSpectrFromSignalSource(int spectrum_size_order)
{
    if (!famplitude_spectrum_computed)
    {
        checkSpectrumSize(static_cast<int>(pow(2.0f, spectrum_size_order)));
        if (!famplitude_spectrum)
            famplitude_spectrum = ippsMalloc_32f(fspectrum_size);
        //QFile f("spectrum.f32");
        //f.open(QIODevice::WriteOnly);
        //f.write((char*)fsource_signal,getSourceSignalSize()*sizeof(float)*2);
        //f.close();
        if (!compute_amp_spectr(getSourceSignal(), getSourceSignalSize(), famplitude_spectrum, spectrum_size_order))
            return nullptr;
        famplitude_spectrum_computed = true;
    }
    return famplitude_spectrum;
}

float* NeuroCalculationEngine::ampSpectrFromSignalFourthDegree(int spectrum_size_order)
{
    if (!famplitude_spectrum_signal_fourth_degree_computed)
    {
        checkSpectrumSize(static_cast<int>(pow(2.0f, spectrum_size_order)));
        if (!famplitude_spectrum_signal_fourth_degree)
            famplitude_spectrum_signal_fourth_degree = ippsMalloc_32f(fspectrum_size);
        if (!compute_amp_spectr(signalFourthDegree(), getSourceSignalSize(), famplitude_spectrum_signal_fourth_degree, spectrum_size_order))
            return nullptr;
        famplitude_spectrum_signal_fourth_degree_computed = true;
    }
    return famplitude_spectrum_signal_fourth_degree;
}

float* NeuroCalculationEngine::ampSpectrFromSignalEighthDegree(int spectrum_size_order)
{
    if (!famplitude_spectrum_signal_eighth_degree_computed)
    {
        checkSpectrumSize(static_cast<int>(pow(2.0f, spectrum_size_order)));
        if (!famplitude_spectrum_signal_eighth_degree)
            famplitude_spectrum_signal_eighth_degree = ippsMalloc_32f(fspectrum_size);
        if (!compute_amp_spectr(signalEighthDegree(), getSourceSignalSize(), famplitude_spectrum_signal_eighth_degree, spectrum_size_order))
            return nullptr;
        famplitude_spectrum_signal_eighth_degree_computed = true;
    }
    return famplitude_spectrum_signal_eighth_degree;
}

float* NeuroCalculationEngine::ampSpectrFromSignalSquare(int spectrum_size_order)
{
    if (!famplitude_spectrum_signal_square_computed)
    {
        checkSpectrumSize(static_cast<int>(pow(2.0f, spectrum_size_order)));
        if (!famplitude_spectrum_signal_square)
            famplitude_spectrum_signal_square = ippsMalloc_32f(fspectrum_size);
        if (!compute_amp_spectr(signalSquare(), getSourceSignalSize(), famplitude_spectrum_signal_square, spectrum_size_order))
            return nullptr;
        famplitude_spectrum_signal_square_computed = true;
    }
    return famplitude_spectrum_signal_square;
}

float* NeuroCalculationEngine::ampSpectrFromEnvelop(int spectrum_size_order)
{
    if (!famplitude_spectrum_envelope_computed)
    {
        checkSpectrumSize(static_cast<int>(pow(2.0f, spectrum_size_order)));
        if (!famplitude_spectrum_envelope)
            famplitude_spectrum_envelope = ippsMalloc_32f(fspectrum_size);
        if (!compute_amp_spectr(envelope(), getSourceSignalSize(), famplitude_spectrum_envelope, spectrum_size_order))
            return nullptr;
        famplitude_spectrum_envelope_computed = true;
    }
    return famplitude_spectrum_envelope;
}


float* NeuroCalculationEngine::ampSpectrFromCorrEnvelope(int spectrum_size_order)
{
    if (!famplitude_spectrum_corr_envelope_computed)
    {
        int sp_size = static_cast<int>(pow(2.0f, spectrum_size_order));
        checkSpectrumSize(sp_size);
        if (!famplitude_spectrum_corr_envelope)
            famplitude_spectrum_corr_envelope = ippsMalloc_32f(fspectrum_size);
        if (!compute_amp_spectr(corrEnvelope(sp_size), sp_size, famplitude_spectrum_corr_envelope, spectrum_size_order))
            return nullptr;
        famplitude_spectrum_corr_envelope_computed = true;
    }
    return famplitude_spectrum_corr_envelope;
}




//!ряд сервисных функций
int     NeuroCalculationEngine::findSignalWideInPoints(float* data, int size, float &max, float &peakSum, float &signalCenter, int max_hole_size, bool putZero)
{
    if (!data)
        return -1;
    //поиск максимума
    int index;
    ippsMaxIndx_32f(data, size, &peakSum, &index);
    max = peakSum;
    int left = index, right = index;
    bool hole_is = max_hole_size;
    max_hole_size = max_hole_size * fspectrum_size / fsample_rate;
    while (left >= 0)
    {
        if (data[left] > max / 3.f)
            peakSum += data[left--];
        else
        {
            if (hole_is)
            {
                int hole = left;
                while ((hole >= 0) && (left - hole < max_hole_size) && (data[hole--] <= max / 3.f))
                {
                }
                if ((hole >= 0) && (data[hole] > max / 3.f) && (left - hole < max_hole_size))
                {
                    hole_is = false;//одной дыры достаточно
                    left = hole;
                    continue;
                }
            }
            break;
        }
    }
    while (right < size)
    {
        if (data[right] > max / 3.f)
            peakSum += data[right++];
        else
        {
            if (hole_is)
            {
                int hole = right;
                while ((hole < size) && (hole - right < max_hole_size) && (data[hole++] <= max / 3.f))
                {
                }
                if ((hole < size) && (data[hole] > max / 3.f) && (hole - right < max_hole_size))
                {
                    hole_is = false;//одной дыры вхатит
                    right = hole;
                    continue;
                }
            }
            break;
        }
    }
    if (putZero)
        ippsZero_32f(&data[left], right - left + 1);
    signalCenter = float(right + left) / 2.0f;
    return (right - left);
}

int     NeuroCalculationEngine::findSignalWideInHz(float* data, int size, float &max, float &peakSum, float &signalCenter, int max_hole_size, bool putZero)
{
    int p = findSignalWideInPoints(data, size, max, peakSum, signalCenter, max_hole_size, putZero);
    if (p < 0)
        return p;
    signalCenter = getFrequency(signalCenter);
    return (p)*fsample_rate / fspectrum_size;
}

//определение ширины сигнала и центра сигнала - возвращается в герцах,  
//max_hole_size - максимальный размер одиночного провала между частями сигнала 
//int       findSignalWideAndZero(float* data, int size,float &max,float &peakSum,float &signalCenter,int max_hole_size=0);

void    NeuroCalculationEngine::normBySumm(float* data, int size)
{
    if (!data)
        return;
    float sum;
    ippsSum_32f(data, size, &sum, ippAlgHintFast);
    ippsDivC_32f_I(sum, data, size);
}

void    NeuroCalculationEngine::meanWindow(float* data, int size, int windowSize)
{
    if ((windowSize < 2) || (windowSize > size / 2))
        return;
    if (!data)
        return;
    //считаем непосредственно
    //ippsSet_32f(1.,data,size);
    float sum;
    ippsSum_32f(data, windowSize - 1, &sum, ippAlgHintFast);
    for (int i = 0; i < size - windowSize; i++)
    {
        sum += data[i + windowSize];
        float temp = sum;
        sum -= data[i];
        data[i] = temp;
    }
    //сдвигаем окно на центр
    ippsMove_32f(data, &data[windowSize / 2], size - windowSize);
    //копируем хвосты окна
    float a = data[windowSize / 2], b = data[size - windowSize / 2];
    for (int i = 0; i < windowSize / 2; i++)
    {
        data[i] = a;
        data[size - i - 1] = b;
    }
    //возвращаем уровень сигнала в исхожное значение
    ippsDivC_32f_I(windowSize, data, size);
}

void    NeuroCalculationEngine::meanWindow(float* source, float *dest, int size, int windowSize)
{
    if ((windowSize < 2) || (windowSize > size / 2))
        return;
    if (!dest)
        return;
    //считаем непосредственно
    //ippsSet_32f(1.,data,size);
    float sum;
    ippsSum_32f(source, windowSize - 1, &sum, ippAlgHintFast);
    float *_dest = &dest[windowSize / 2];
    for (int i = 0; i < size - windowSize; i++)
    {
        sum += source[i + windowSize];
        float temp = sum;
        sum -= source[i];
        _dest[i] = temp;
    }
    //копируем хвосты окна
    float a = dest[windowSize / 2], b = dest[size - windowSize / 2];
    for (int i = 0; i < windowSize / 2; i++)
    {
        dest[i] = a;
        dest[size - i - 1] = b;
    }
    //возвращаем уровень сигнала в исхожное значение
    ippsDivC_32f_I(windowSize, dest, size);
}


int    NeuroCalculationEngine::getSpectrumOrder(void)
{
    return getSpectrumOrder(fspectrum_size);
};

int    NeuroCalculationEngine::getSpectrumOrder(int spectrum_size)
{
    int order = 0;
    while (spectrum_size > 1)
    {
        spectrum_size /= 2;
        order++;
    }
    return order;
}


//получение номера спектрального коэффициента от частоты 
int     NeuroCalculationEngine::getSpectrPos(float frequency)
{
    if (frequency <= -fsample_rate / 2.f)
    {
        Q_ASSERT(false);
        frequency = -fsample_rate / 2;
    }
    if (frequency >= fsample_rate / 2.f)
    {
        Q_ASSERT(false);
        frequency = fsample_rate / 2;
    }
    return static_cast<int>((frequency + fsample_rate / 2.f) / fsample_rate * fspectrum_size + 0.5f);
}

//получение частоты от спектрального коэффициента
float   NeuroCalculationEngine::getFrequency(float spectrPos)
{
    if (spectrPos < 0.f)
    {
        Q_ASSERT(false);
        spectrPos = 0.0;
    }
    if (spectrPos >= fspectrum_size)
    {
        Q_ASSERT(false);
        spectrPos = fspectrum_size - 1;
    }
    return (spectrPos)*fsample_rate / fspectrum_size - fsample_rate / 2.f;
}

int    NeuroCalculationEngine::getSpectrumSize(void)
{
    return fspectrum_size;
};

int    NeuroCalculationEngine::getSpectrumSize(int spectrum_size_order)
{
    return static_cast<int>(pow(2.0f, spectrum_size_order) + 0.5);
}


float* NeuroCalculationEngine::freqDetector(int offset)
{
    if (!ffreq_detector_computed)
    {
        if (!ffreq_detector)
            ffreq_detector = ippsMalloc_32f(fsource_max_signal_size);
        //float *temp = tempArray();
        const Ipp32fc *signal = fsource_signal;//getSourceSignal();
        int signal_size = fsource_signal_size;//getSourceSignalSize();
//              DEM[i]=     fais_temp_buffer[i+step].re*fais_temp_buffer[i].im-
//              fais_temp_buffer[i+step].im*fais_temp_buffer[i].re;
        /*
        float *re = temp;
        float *im = &temp[signal_size];
        ippsCplxToReal_32fc(signal,re,im,signal_size);
        ippsMul_32f(im,&re[offset],&ffreq_detector[ftail_signal_size],signal_size-offset);//im*re(+offset)
        ippsMul_32f_I(&im[offset],re,signal_size-offset);//re*im(+offset)
        ippsSub_32f_I(re,&ffreq_detector[ftail_signal_size],signal_size-offset);
        */
        for (int i = 0; i < signal_size - offset; i++)
            ffreq_detector[i + ftail_signal_size] = signal[i + offset].re*signal[i].im -
            signal[i + offset].im*signal[i].re;
        ffreq_detector_computed = true;
    }
    return ffreq_detector;
}

//!более усложненный вариант частотного детектора - пригоден для анализа многочастотных сигналов
float* NeuroCalculationEngine::freqDetector2(void)
{
    if (!ffreq_detector_computed)
    {
        if (!ffreq_detector)
            ffreq_detector = ippsMalloc_32f(fsource_max_signal_size);
        ffreq_detector_computed = true;
        const Ipp32fc *signal = /*fsource_signal;//*/getSourceSignal();
        int signal_size = /*fsource_signal_size;*/getSourceSignalSize();
        Ipp32fc* temp = reinterpret_cast<Ipp32fc*>(tempArray());
        ippsDiv_32fc(&signal[0], &signal[1], &temp[1], signal_size - 1);
        if ((fabs(previos_signal_point.re) < 0.000000001f) && (fabs(previos_signal_point.im) < 0.00000001f))
            temp[0] = temp[1];//чтобы избежать деление на ноль
        else
        {
            float t = signal[0].re*signal[0].re + signal[0].im*signal[0].im;
            temp[0].re = (previos_signal_point.re*signal[0].re - previos_signal_point.im*signal[0].im) / t;
            temp[0].im = (previos_signal_point.re*signal[0].im - previos_signal_point.im*signal[0].re) / t;
        }
        previos_signal_point = signal[signal_size - 1];
        ippsPhase_32fc(temp, &ffreq_detector[/*ftail_signal_size*/0], signal_size);
        ippsMulC_32f_I(static_cast<float>(fsample_rate / 2. / IPP_PI), &ffreq_detector[/*ftail_signal_size*/0], signal_size);
    }
    return ffreq_detector;
}
/*
//!третий вариант частотного детектора - предложенный Казаковым Антоном - надо тестить
float* CalculationEngine::freqDetector3(void)
{
    if(!ffreq_detector_computed)
    {   if(!ffreq_detector)
            ffreq_detector = ippsMalloc_32f(fsource_max_signal_size);
        ffreq_detector_computed = true;
        const Ipp32fc *signal = getSourceSignal();
        int signal_size = getSourceSignalSize();
        Ipp32fc* temp = (Ipp32fc*)tempArray();
        ippsDiv_32fc(&signal[0],&signal[1],&temp[1],signal_size-1);
        if((fabs(previos_signal_point.re)<0.000000001)&&(fabs(previos_signal_point.im)<0.00000001))
            temp[0] = temp[1];//чтобы избежать деление на ноль
        else
        {   float t = signal[0].re*signal[0].re+signal[0].im*signal[0].im;
            temp[0].re = (previos_signal_point.re*signal[0].re-previos_signal_point.im*signal[0].im)/t;
            temp[0].im = (previos_signal_point.re*signal[0].im-previos_signal_point.im*signal[0].re)/t;
        }
        previos_signal_point = signal[signal_size-1];
        ippsPhase_32fc(temp,&ffreq_detector[0],signal_size);
        ippsMulC_32f_I(fsample_rate/2./IPP_PI,&ffreq_detector[0],signal_size);
    }
    return ffreq_detector;
}
*/

//!возвращает исходный радиосигнал
const Ipp32fc   *NeuroCalculationEngine::getSourceSignal(void)
{
    return fsource_signal;
}

//!получение размера хвоста от предыдущего выполнения процедуры append
int NeuroCalculationEngine::getTailSignalSize(void)
{
    return ftail_signal_size;
};

int     NeuroCalculationEngine::getSampleRate(void)
{
    return fsample_rate;
};

int     NeuroCalculationEngine::getBand(void)
{
    return fband;
};

int     NeuroCalculationEngine::getSourceSignalSize(void)
{
    return fsource_signal_size + ftail_signal_size;
};

void NeuroCalculationEngine::signalToFloat(const Ipp32fc* fSignal, float* resultSignal, int size_of_signal)
{
    for (int i = 0, j = 0; i < size_of_signal * 2, j < size_of_signal; i = i + 2, j = j++)
    {
        resultSignal[i] = fSignal[j].re;
        resultSignal[i + 1] = fSignal[j].im;
    }
}

void NeuroCalculationEngine::invertSp(float* nSpectr, int SizeOfSp)
{
    int size = 0;
    float temp_float = 0.;

    for (unsigned int i = 0, size = (SizeOfSp >> 1); size > i; i++)
    {
        temp_float = nSpectr[size + i];
        nSpectr[size + i] = nSpectr[i];
        nSpectr[i] = temp_float;
    }
}

Ipp32fc* NeuroCalculationEngine::compute_auto_corr_tfc(Ipp32fc* signal, Ipp32fc* autocorr, int power)
{
    ippsZero_32fc(ftemp_sig_array_tfc, fsource_signal_size);
    Ipp32fc* tfft_arr;
    tfft_arr = fftArr(power);

    for (int i = 0; i < fsource_signal_size; i++)
    {
        ftemp_sig_array_tfc[i].re = sqrt((tfft_arr[i].re * tfft_arr[i].re) + (tfft_arr[i].im * tfft_arr[i].im));
        ftemp_sig_array_tfc[i].im = 0.0;
    }

    if (power != _fft_old->Order())
        _fft_old->Order(power);

    _fft_old->RunInversion(tfft_arr, autocorr);
    fautocorr_signal_source_tfc_computed = true;

    return autocorr;
}

Ipp32fc* NeuroCalculationEngine::compute_auto_corr_tfc_sec(Ipp32fc* signal, Ipp32fc* autocorr, int power)
{
    /*ippsZero_32fc(ftemp_sig_array_tfc, fsource_signal_size);
    Ipp32fc* tfft_arr;
    tfft_arr = fftArr(power);*/
    //int shift = 40; //константа сдвига массива для проведения корреляции
    float sumIm = 0.0;
    float sumRe = 0.0;
    //float sum = 0.0;

    for (int i = 0; i < fsource_signal_size; i++)
    {
        sumIm = 0;
        sumRe = 0;
        //sum = 0;
        for (int j = 0; j < fsource_signal_size - i; j++)
        {
            sumRe += signal[j].re * signal[j + i].re;
            sumIm += signal[j].im * signal[j + i].im;
            //sumRe += (In[j].re * In[j+i].re) - (In[j].im * In[j+i].im);
            //sumIm += (In[j].re * In[j+i].im) + (In[j+i].re * In[j].im);
            //sum += ((In[j].re  * In[j+i].re) - (In[j].im * In[j+i].im)) + ((In[j].re * In[j+i].im) + (In[j+i].re * In[j].im));

        }
        autocorr[i].re = sumRe / (fsource_signal_size - i);
        autocorr[i].im = sumIm / (fsource_signal_size - i);
    }



    fautocorr_signal_source_tfc_computed = true;

    return autocorr;
}

int NeuroCalculationEngine::findAnotherPeaksWithSpace(float* array, int startpoint, int endpoint, int sec_startpoint, int sec_endpoint, float summ, float k)
{
    int i = 0, max = 0;
    int num_max = 0;
    for (i = startpoint, max = 0; i < endpoint; i++)
    {
        if ((array[i] > summ *k) && (i<(sec_startpoint) || i>(sec_endpoint)))
        {
            summ = array[i];
            max = i;
            num_max++;
        } /* max */
    }
    return num_max;
}

float NeuroCalculationEngine::findEnergy(const Ipp32fc* sig, int min_spectrum_size, int startpoint)
{
    float summ = 0, energy = 0;

    for (int i = startpoint; i < startpoint + min_spectrum_size; i++)
    {
        summ = summ + (sig[i].re * sig[i].re + sig[i].im * sig[i].im);
    }
    energy = summ / min_spectrum_size;
    return energy;
}


bool NeuroCalculationEngine::energyTest(const Ipp32fc* sig, int miniSize, float minEnergy, float maxEnergy, int startpoint)
{
    if (findEnergy(sig, miniSize, startpoint) > (minEnergy + (maxEnergy - minEnergy) / 5))
        return true;
    else
        return false;
}

void NeuroCalculationEngine::findOffsetAndBand(float* spectrum, int &band, int &offset, int spectrum_size_order)
{
    int leftgr = 0;
    int rightgr = 0;
    int kpor = 150; // 0,15 от уровня сигнала
    float max = 0.;
    float hp = 0.;
    int maxofarray;
    //int SizeofSp = pow(2,SpectrumSizeOrder);
    int _size_of_spectrum = pow(2, spectrum_size_order);
    ippsZero_32f(ftemp_amplitude_spectrum, _size_of_spectrum);
    ippsCopy_32f(spectrum, ftemp_amplitude_spectrum, _size_of_spectrum);
    for (int i = 10; i < _size_of_spectrum - 10; i++)
    {
        if (ftemp_amplitude_spectrum[i] > 0.)
        {
            ftemp_amplitude_spectrum[i] = (ftemp_amplitude_spectrum[i] + ftemp_amplitude_spectrum[i + 1] + ftemp_amplitude_spectrum[i + 2] + ftemp_amplitude_spectrum[i + 3])*0.25;
            if (ftemp_amplitude_spectrum[i] > max)
            {
                max = ftemp_amplitude_spectrum[i];
                maxofarray = i;
            }
        }
    }
    if (max > 0.)
    {
        max = 1000. / max;
    }
    for (int i = 10; i < _size_of_spectrum - 10; i++)
    {
        ftemp_amplitude_spectrum[i] = ftemp_amplitude_spectrum[i] * max;
    }

    for (int i = 10; i < _size_of_spectrum - 10; i++)
    {
        if (ftemp_amplitude_spectrum[i] > kpor)
        {
            leftgr = i;
            break;
        }
    }
    for (int i = _size_of_spectrum - 10; i > 10; i--)
    {
        if (ftemp_amplitude_spectrum[i] > kpor)
        {
            rightgr = i + 1;
            break;
        }
    }
    band = abs(rightgr - leftgr);
    band = band * (fsample_rate / _size_of_spectrum);
    hp = (float)(leftgr + rightgr);
    hp = hp * 0.5;
    offset = (int)((hp - (_size_of_spectrum / 2)) * (float)fsample_rate / (float)_size_of_spectrum);
}

int NeuroCalculationEngine::findMaxInRegion(float* Array, int startpoint, int endpoint, float summ, float k)
{
    int i;
    double maxofSp;
    int max;
    for (i = startpoint, max = 0, maxofSp = summ * k; i < endpoint; i++)
    {
        if (Array[i] > maxofSp)
        {
            maxofSp = Array[i];
            max = i;
        }
    }
    return max;
}

Ipp32fc* NeuroCalculationEngine::autoCorr(int spectrum_size_order)
{
    if (0 == fsource_signal)
        return 0;
    if (!fautocorr_signal_source_tfc_computed)
    {
        /*if (0 == compute_auto_corr_tfc(fsource_signal, fautocorr_signal_source_tfc, spectrum_size_order))*/
        if (0 == compute_auto_corr_tfc_sec(fsource_signal, fautocorr_signal_source_tfc, spectrum_size_order))
            return 0;
        fautocorr_signal_source_tfc_computed = true;
        //теперь заодно и формируем корреляционную функцию формата kls
        signalToFloat(fautocorr_signal_source_tfc, fautocorr_klsformat, fsource_signal_size);
        ippsAbs_32f(fautocorr_klsformat, fautocorr_klsformat, fsource_signal_size);
        fautocorr_klsformat_computed = true;
    }

    return fautocorr_signal_source_tfc;
}

float* NeuroCalculationEngine::autoCorrKlsFormat()
{
    if (!fautocorr_klsformat_computed)
        return 0;
    else
        return fautocorr_klsformat;
}

Ipp32fc* NeuroCalculationEngine::fftArr(int spectrum_size_order)
{
    if (0 == fsource_signal)
        return 0;
    if (!ffft_arr_computed)
    {
        if (0 == compute_fft_arr(fsource_signal, ffft_arr, spectrum_size_order))
            return 0;
        ffft_arr_computed = true;
    }
    return ffft_arr;
}

Ipp32fc* NeuroCalculationEngine::compute_fft_arr(Ipp32fc* signal, Ipp32fc* fftarr, int power)
{
    ippsZero_32fc(fftarr, fsource_signal_size);
    if (power != _fft_old->Order())
        _fft_old->Order(power);

    _fft_old->Run((TFloatComplex*)signal, fftarr);
    return fftarr;
}

float* NeuroCalculationEngine::ampSpectrFromAmpSpectrFromSignalSource(int spectrum_size_order)
{
    int _size_of_spectrum = pow(2, spectrum_size_order);
    if (!famplitude_spectrum_amplitude_spectrum_signal_source_computed)
    {
        checkSpectrumSize(static_cast<int>(pow(2.0f, spectrum_size_order)));
        if(!famplitude_spectrum_amplitude_spectrum_signal_source)
            famplitude_spectrum_amplitude_spectrum_signal_source = ippsMalloc_32f(fspectrum_size);
        if (!compute_amp_spectr(ampSpectrFromSignalSource(spectrum_size_order), _size_of_spectrum,famplitude_spectrum_amplitude_spectrum_signal_source, spectrum_size_order))
            return nullptr;
        famplitude_spectrum_amplitude_spectrum_signal_source_computed = true;
    }
    return famplitude_spectrum_amplitude_spectrum_signal_source;
}

float* NeuroCalculationEngine::ampSpectrFromAmpSpectrFromSignalSquare(int spectrum_size_order)
{
    int _size_of_spectrum = pow(2, spectrum_size_order);
    if (!famplitude_spectrum_amplitude_spectrum_signal_square_computed)
    {
        checkSpectrumSize(static_cast<int>(pow(2.0f, spectrum_size_order)));
        if (!famplitude_spectrum_amplitude_spectrum_signal_square)
            famplitude_spectrum_amplitude_spectrum_signal_square = ippsMalloc_32f(fspectrum_size);
        if (!compute_amp_spectr(ampSpectrFromSignalSquare(spectrum_size_order), _size_of_spectrum, famplitude_spectrum_amplitude_spectrum_signal_square, spectrum_size_order))
            return nullptr;
        famplitude_spectrum_amplitude_spectrum_signal_square_computed = true;
    }
    return famplitude_spectrum_amplitude_spectrum_signal_square;
}

float* NeuroCalculationEngine::ampSpectrFromAmpSpectrFromSignalFourthDegree(int spectrum_size_order)
{
    int _size_of_spectrum = pow(2, spectrum_size_order);
    if (!famplitude_spectrum_amplitude_spectrum_signal_fourth_degree_computed)
    {
        checkSpectrumSize(static_cast<int>(pow(2.0f, spectrum_size_order)));
        if (!famplitude_spectrum_amplitude_spectrum_signal_fourth_degree)
            famplitude_spectrum_amplitude_spectrum_signal_fourth_degree = ippsMalloc_32f(fspectrum_size);
        if (!compute_amp_spectr(ampSpectrFromSignalFourthDegree(spectrum_size_order), _size_of_spectrum, famplitude_spectrum_amplitude_spectrum_signal_fourth_degree, spectrum_size_order))
            return nullptr;
        famplitude_spectrum_amplitude_spectrum_signal_fourth_degree_computed = true;
    }
    return famplitude_spectrum_amplitude_spectrum_signal_fourth_degree;
}

float* NeuroCalculationEngine::ampSpectrFromAmpSpectrFromSignalEightDegree(int spectrum_size_order)
{
    int _size_of_spectrum = pow(2, spectrum_size_order);
    if (!famplitude_spectrum_amplitude_spectrum_signal_eight_degree_computed)
    {
        checkSpectrumSize(static_cast<int>(pow(2.0f, spectrum_size_order)));
        if (!famplitude_spectrum_amplitude_spectrum_signal_eight_degree)
            famplitude_spectrum_amplitude_spectrum_signal_eight_degree = ippsMalloc_32f(fspectrum_size);
        if (!compute_amp_spectr(ampSpectrFromSignalEighthDegree(spectrum_size_order), _size_of_spectrum, famplitude_spectrum_amplitude_spectrum_signal_eight_degree, spectrum_size_order))
            return nullptr;
        famplitude_spectrum_amplitude_spectrum_signal_eight_degree_computed = true;
    }
    return famplitude_spectrum_amplitude_spectrum_signal_eight_degree;
}

float* NeuroCalculationEngine::ampSpectrFromEnvelopInv(int spectrum_size_order)
{
    if (!famplitude_spectrum_envelope_inv_computed)
    {
        checkSpectrumSize(static_cast<int>(pow(2.0f, spectrum_size_order)));
        if (!famplitude_spectrum_envelope_inv)
            famplitude_spectrum_envelope_inv = ippsMalloc_32f(fspectrum_size);
        if (!compute_amp_spectr(envelope(), getSourceSignalSize(), famplitude_spectrum_envelope_inv, spectrum_size_order))
            return nullptr;
        invertSp(famplitude_spectrum_envelope_inv, fspectrum_size);
        famplitude_spectrum_envelope_inv_computed = true;
    }

    return famplitude_spectrum_envelope_inv;
}


float* NeuroCalculationEngine::ampSpectrFromModDerEnvelope(int spectrum_size_order)
{
    if (!famplitude_spectrum_mod_der_envelope_computed)
    {
        checkSpectrumSize(static_cast<int>(pow(2.0f, spectrum_size_order)));
        if (!famplitude_spectrum_mod_der_envelope)
            famplitude_spectrum_mod_der_envelope = ippsMalloc_32f(fspectrum_size);
        if (!compute_amp_spectr(modDerEnvelope(), getSourceSignalSize(), famplitude_spectrum_mod_der_envelope, spectrum_size_order))
            return nullptr;
        famplitude_spectrum_mod_der_envelope_computed = true;
    }

    return famplitude_spectrum_mod_der_envelope;
}



float*	NeuroCalculationEngine::compute_amp_spectr_from_amp_spectr(const float* signal, int signal_size, float *spectrum, int power)
{
    int spectrum_size = static_cast<int>(pow(2.0f, power));
    if (fspectrum_size > signal_size)
        return nullptr;
    if (!ftemp_amplitude_spectrum)
        ftemp_amplitude_spectrum = ippsMalloc_32f(fspectrum_size);
    _r_fft->magnutude(signal, spectrum_size, spectrum);
    for (int i = fspectrum_size; i < signal_size - fspectrum_size; i += fspectrum_size)
    {
        _r_fft->magnutude(&signal[i], spectrum_size, ftemp_amplitude_spectrum);
        ippsAdd_32f_I(ftemp_amplitude_spectrum, spectrum, fspectrum_size);
    }
    normBySumm(spectrum, fspectrum_size);
    return spectrum;
 

    
   
}

float* NeuroCalculationEngine::modDerEnvelope(void)
{
    if (!fsource_signal)
        return nullptr;
    if (!fmod_der_envelope_computed)
    {
        if (!fmod_der_envelope)
            fmod_der_envelope = ippsMalloc_32f(fsource_max_signal_size);
        if (!fenvelope_computed)
        {
            envelope();
        }
        for (int i = 0; i < fsource_signal_size - 1; i++)
        {
            fmod_der_envelope[i] = fenvelope[i] - fenvelope[i + 1];
        }
        ippsAbs_32f(fmod_der_envelope, fmod_der_envelope, fsource_signal_size);
        fmod_der_envelope_computed = true;
    }
    return fmod_der_envelope;
}

float* NeuroCalculationEngine::powerSpectrum2avModDerFreq(int spectrum_size_order)
{
    if (!fpower_spectrum_2av_mod_der_freq_computed)
    {
        checkSpectrumSize(static_cast<int>(pow(2.0f, spectrum_size_order)));
        if (!fpower_spectrum_2av_mod_der_freq)
            fpower_spectrum_2av_mod_der_freq = ippsMalloc_32f(fspectrum_size);
        if (!compute_amp_spectr(avModDerFreq(), getSourceSignalSize(), fpower_spectrum_2av_mod_der_freq, spectrum_size_order))
            return nullptr;
        fpower_spectrum_2av_mod_der_freq_computed = true;
    }
    return fpower_spectrum_2av_mod_der_freq;
}

float* NeuroCalculationEngine::powerSpectrum10avFreq(int spectrum_size_order)
{
    if (!fpower_spectrum_10av_freq_computed)
    {
        checkSpectrumSize(static_cast<int>(pow(2.0f, spectrum_size_order)));
        if (!fpower_spectrum_10av_freq)
            fpower_spectrum_10av_freq = ippsMalloc_32f(fspectrum_size);
        if (!compute_amp_spectr(avFreq(), getSourceSignalSize(), fpower_spectrum_10av_freq, spectrum_size_order))
            return nullptr;
        fpower_spectrum_10av_freq_computed = true;
    }
    return famplitude_spectrum_envelope;
}

float* NeuroCalculationEngine::avFreq(void)
{
    if (!fsource_signal)
        return nullptr;
    if (!av_freq_computed)
    {
        if (!av_freq)
            av_freq = ippsMalloc_32f(fsource_max_signal_size);
        if (!ffrequancy_computed)
        {
            frequency();
        }
        ippsCopy_32f(ffrequancy, av_freq, fsource_max_signal_size);
        av_freq = computeAverage(av_freq, fsource_max_signal_size, 10);
        av_freq_computed = true;
    }
    return av_freq;
}

float* NeuroCalculationEngine::avModDerFreq(void)
{
    if (!fsource_signal)
        return nullptr;
    if (!av_mod_der_freq_computed)
    {
        if (!av_mod_der_freq)
            av_mod_der_freq = ippsMalloc_32f(fsource_max_signal_size);
        if (!ffrequancy_computed)
        {
            frequency();
        }
        for (int i = 0; i < fsource_signal_size - 1; i++)
        {
            av_mod_der_freq[i] = ffrequancy[i] - ffrequancy[i + 1];
        }
        ippsAbs_32f(av_mod_der_freq, av_mod_der_freq, fsource_signal_size);
        av_mod_der_freq = computeAverage(av_mod_der_freq, fsource_max_signal_size, 2);
        av_mod_der_freq_computed = true;
    }
    return av_mod_der_freq;
}


float* NeuroCalculationEngine::computeAverage(float* signal, int signal_size, int average_number)
{
    for (int i = 0; i < signal_size - average_number; i++)
    {
        for (int j = 0; j < average_number; j++)
        {
            signal[i] += signal[i + j];
        }
        signal[i] = signal[i] / average_number;
    }
    return signal;
}

float* NeuroCalculationEngine::normAmpSpectrum(float* spectrum, int coeff, int spectrum_size_order)
{
    int spectrum_size = static_cast<int>(pow(2.0f, spectrum_size_order));
    float maxOfSp = 0.;
    int maxInd = 0;
    float tmp_val = 0.;
    if (!fsource_signal)
        return nullptr;
    if (!famplitude_norm_spectrum_computed)
    {
        if (!famplitude_norm_spectrum)
            famplitude_norm_spectrum = ippsMalloc_32f(fsource_max_signal_size);

        ippsMaxIndx_32f(reinterpret_cast<const float*>(spectrum), spectrum_size, &maxOfSp, &maxInd);
        if (maxOfSp > 0.)
        {
            tmp_val = coeff / maxOfSp;
        }
        for (int i = 0; i < spectrum_size; i++)
        {
            famplitude_norm_spectrum[i] = spectrum[i] * tmp_val;
        }
        ippsMulC_32f(spectrum, tmp_val, famplitude_norm_spectrum, spectrum_size);
        famplitude_norm_spectrum_computed = true;
    }
    return famplitude_norm_spectrum;
}

