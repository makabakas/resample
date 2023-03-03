clear
clc

num_zeros_ = 10;
tar_fs = 8000;

[y, fs] = audioread('C:\Users\tao.li\Downloads\1221-135766-0002.wav');
y = y';

% cutoff frequency is min(orig_fs, tar_fs), because fh is 0.5*fs.
filter_cutoff_ = 0.5 * 0.99 * min(tar_fs, fs); 
base_freq = gcd(fs, tar_fs);

input_samples_in_unit_ = fs / base_freq;
output_samples_in_unit_ = tar_fs / base_freq;

%% set Indexes and calculate filter
first_index_ = zeros(output_samples_in_unit_, 1);
weights_ = cell(output_samples_in_unit_, 1);
window_width = num_zeros_ / (2 * filter_cutoff_);

for i = 0:output_samples_in_unit_-1
    output_t = i / tar_fs;
    min_t = output_t - window_width;
    max_t = output_t + window_width;
    
    min_input_index = ceil(min_t * fs);
    max_input_index = floor(max_t * fs);
    
    num_indices = max_input_index - min_input_index + 1;
    
    first_index_(i+1) = min_input_index;
    
    weight = zeros(num_indices, 1);
    for j = 0: num_indices-1
        input_index = min_input_index + j;
        input_t = input_index / fs;
        delta_t = input_t - output_t;
        
        if(abs(delta_t) < num_zeros_ / (2.0 * filter_cutoff_))
            window = 0.5 * (1 + cos(2 * pi * filter_cutoff_ / num_zeros_ * delta_t));
        else
            window = 0;
        end
        
        if(delta_t ~= 0)
            filter = sin(2 * pi * filter_cutoff_ * delta_t) / (pi * delta_t);
        else
            filter = 2 * filter_cutoff_;
        end
        
        weight(j+1) = window * filter / fs;
    end
    weights_{i+1} = weight;
end
input_sample_offset_ = 0;
output_sample_offset_ = 0;
input_remainder_ = [];

input_dim = length(y);
flush = false;

tot_input_samp = input_sample_offset_ + input_dim;

tick_freq = lcm(fs, tar_fs);
ticks_per_input_period = tick_freq / fs;

interval_length_in_ticks = tot_input_samp * ticks_per_input_period;
if ~flush
    window_width = num_zeros_ / (2.0 * filter_cutoff_);
    window_width_ticks = floor(window_width * tick_freq);
    interval_length_in_ticks = interval_length_in_ticks - window_width_ticks;
end

ticks_per_output_period = tick_freq / tar_fs;
last_output_samp = interval_length_in_ticks / ticks_per_output_period;
if (last_output_samp * ticks_per_output_period == interval_length_in_ticks)
    last_output_samp = last_output_samp - 1;
end
tot_output_samp = last_output_samp + 1;

output = zeros(tot_output_samp - output_sample_offset_, 1);

for samp_out = output_sample_offset_ : tot_output_samp - 1
    unit_index = idivide(int32(samp_out), int32(output_samples_in_unit_));
    samp_out_wrapped = int32(samp_out - unit_index * output_samples_in_unit_);
    first_samp_in = first_index_(samp_out_wrapped + 1) + unit_index * input_samples_in_unit_;
    
    weights = weights_{samp_out_wrapped + 1};
    
    first_input_index = first_samp_in - input_sample_offset_;
    
    if (first_input_index >= 0) && (first_input_index + length(weights) <= input_dim)
        this_output = y(first_input_index+1 : first_input_index + length(weights)) * weights;
    else
        this_output = 0;
        for i = 0 : length(weights)-1
            weight = weights(i+1);
            input_index = first_input_index + i;
            if(input_index < 0) && (length(input_remainder_) + input_index >= 0)
                this_output = this_output + ...
                        weight * input_remainder_(length(input_remainder_) ... 
                        + input_index);
                    
            elseif (input_index >= 0) && (input_index < input_dim)
                this_output = this_output + weight * y(input_index+1);
            end
        end
    end
    output_index = samp_out - output_sample_offset_ + 1;
    output(output_index) = this_output;
end

audiowrite('1221-135766-0002.wav', output, tar_fs);