#include "kseq/kseq.h"
#include "common.h"

// Optimized CUDA kernel for string matching
__global__ void matchStringsKernel(
    const char *__restrict__ d_sample_sequences,
    const int *__restrict__ d_sample_offsets,
    const int *__restrict__ d_sample_lengths,
    const char *__restrict__ d_sample_qualities,
    const char *__restrict__ d_signature_sequences,
    const int *__restrict__ d_signature_offsets,
    const int *__restrict__ d_signature_lengths,
    const int num_samples,
    const int num_signatures,
    double *__restrict__ d_match_matrix)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total_pairs = num_samples * num_signatures;

    if (idx < total_pairs)
    {
        int sample_idx = idx / num_signatures;
        int signature_idx = idx % num_signatures;

        int sample_start = d_sample_offsets[sample_idx];
        const char *sample_seq = &d_sample_sequences[sample_start];
        const char *sample_qual = &d_sample_qualities[sample_start];
        int sample_len = d_sample_lengths[sample_idx];

        int signature_start = d_signature_offsets[signature_idx];
        const char *signature_seq = &d_signature_sequences[signature_start];
        int signature_len = d_signature_lengths[signature_idx];

        double match_score = -1.0;

        for (int i = 0; i <= sample_len - signature_len; i++)
        {
            bool match_found = true;
            double match_score_local = 0.0;

            for (int j = 0; j < signature_len; j++)
            {
                char sample_char = sample_seq[i + j];
                char signature_char = signature_seq[j];
                char sample_qual_char = sample_qual[i + j];

                if (sample_char != signature_char &&
                    sample_char != 'N' &&
                    signature_char != 'N')
                {
                    match_found = false;
                    break;
                }
                match_score_local += (sample_qual_char - 33);
            }
            if (match_found)
            {
                match_score = match_score_local / signature_len;
                break;
            }
        }

        d_match_matrix[sample_idx * num_signatures + signature_idx] = match_score;
    }
}

void runMatcher(const std::vector<klibpp::KSeq> &samples,
                const std::vector<klibpp::KSeq> &signatures,
                std::vector<MatchResult> &matches)
{
    int num_samples = samples.size();
    int num_signatures = signatures.size();

    // Allocate host arrays for lengths and offsets
    int *h_sample_lengths = new int[num_samples];
    int *h_signature_lengths = new int[num_signatures];
    int *h_sample_offsets = new int[num_samples];
    int *h_signature_offsets = new int[num_signatures];

    // Calculate total lengths and offsets
    size_t total_sample_length = 0;
    for (int i = 0; i < num_samples; ++i)
    {
        h_sample_lengths[i] = samples[i].seq.size();
        h_sample_offsets[i] = total_sample_length;
        total_sample_length += h_sample_lengths[i];
    }

    size_t total_signature_length = 0;
    for (int i = 0; i < num_signatures; ++i)
    {
        h_signature_lengths[i] = signatures[i].seq.size();
        h_signature_offsets[i] = total_signature_length;
        total_signature_length += h_signature_lengths[i];
    }

    // Allocate concatenated sequences and qualities
    char *h_sample_sequences = new char[total_sample_length];
    char *h_sample_qualities = new char[total_sample_length];
    char *h_signature_sequences = new char[total_signature_length];

    // Copy sequences and qualities into concatenated arrays
    for (int i = 0; i < num_samples; ++i)
    {
        memcpy(&h_sample_sequences[h_sample_offsets[i]],
               samples[i].seq.c_str(), h_sample_lengths[i]);
        memcpy(&h_sample_qualities[h_sample_offsets[i]],
               samples[i].qual.c_str(), h_sample_lengths[i]);
    }

    for (int i = 0; i < num_signatures; ++i)
    {
        memcpy(&h_signature_sequences[h_signature_offsets[i]],
               signatures[i].seq.c_str(), h_signature_lengths[i]);
    }

    // Allocate device memory
    char *d_sample_sequences, *d_sample_qualities;
    int *d_sample_offsets, *d_sample_lengths;
    char *d_signature_sequences;
    int *d_signature_offsets, *d_signature_lengths;
    double *d_match_matrix;

    cudaMalloc(&d_sample_sequences, total_sample_length * sizeof(char));
    cudaMalloc(&d_sample_qualities, total_sample_length * sizeof(char));
    cudaMalloc(&d_sample_offsets, num_samples * sizeof(int));
    cudaMalloc(&d_sample_lengths, num_samples * sizeof(int));

    cudaMalloc(&d_signature_sequences, total_signature_length * sizeof(char));
    cudaMalloc(&d_signature_offsets, num_signatures * sizeof(int));
    cudaMalloc(&d_signature_lengths, num_signatures * sizeof(int));

    cudaMalloc(&d_match_matrix, num_samples * num_signatures * sizeof(double));

    // Copy data to device
    cudaMemcpy(d_sample_sequences, h_sample_sequences,
               total_sample_length * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sample_qualities, h_sample_qualities,
               total_sample_length * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sample_offsets, h_sample_offsets,
               num_samples * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sample_lengths, h_sample_lengths,
               num_samples * sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(d_signature_sequences, h_signature_sequences,
               total_signature_length * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_signature_offsets, h_signature_offsets,
               num_signatures * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_signature_lengths, h_signature_lengths,
               num_signatures * sizeof(int), cudaMemcpyHostToDevice);

    // Configure kernel launch parameters
    int total_pairs = num_samples * num_signatures;
    int threadsPerBlock = 256;
    int blocksPerGrid = (total_pairs + threadsPerBlock - 1) / threadsPerBlock;

    // Launch the optimized kernel
    matchStringsKernel<<<blocksPerGrid, threadsPerBlock>>>(
        d_sample_sequences,
        d_sample_offsets,
        d_sample_lengths,
        d_sample_qualities,
        d_signature_sequences,
        d_signature_offsets,
        d_signature_lengths,
        num_samples,
        num_signatures,
        d_match_matrix);

    // Copy match matrix back to host
    double *h_match_matrix = new double[num_samples * num_signatures];
    cudaMemcpy(h_match_matrix, d_match_matrix,
               num_samples * num_signatures * sizeof(double), cudaMemcpyDeviceToHost);

    // Process the match results
    for (int i = 0; i < num_samples; ++i)
    {
        for (int j = 0; j < num_signatures; ++j)
        {
            if (h_match_matrix[i * num_signatures + j] != -1)
            {
                MatchResult result;
                result.sample_name = samples[i].name;
                result.signature_name = signatures[j].name;
                result.match_score = h_match_matrix[i * num_signatures + j];
                matches.push_back(result);
            }
        }
    }

    // Clean up host and device memory
    delete[] h_sample_lengths;
    delete[] h_signature_lengths;
    delete[] h_sample_offsets;
    delete[] h_signature_offsets;
    delete[] h_sample_sequences;
    delete[] h_sample_qualities;
    delete[] h_signature_sequences;
    delete[] h_match_matrix;

    cudaFree(d_sample_sequences);
    cudaFree(d_sample_qualities);
    cudaFree(d_sample_offsets);
    cudaFree(d_sample_lengths);
    cudaFree(d_signature_sequences);
    cudaFree(d_signature_offsets);
    cudaFree(d_signature_lengths);
    cudaFree(d_match_matrix);
}