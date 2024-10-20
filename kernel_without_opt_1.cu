#include "kseq/kseq.h"
#include "common.h"

__global__ void matchStringsKernel(
    char **__restrict__ d_sample_sequences,
    char **__restrict__ d_sample_qualities,
    const int *__restrict__ d_sample_lengths,
    char **__restrict__ d_signature_sequences,
    const int *__restrict__ d_signature_lengths,
    const int num_samples,
    const int num_signatures,
    double *__restrict__ d_match_matrix)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total_pairs = num_samples * num_signatures;

    if (idx < total_pairs)
    {
        // Get the sample and signature indices
        int sample_idx = idx / num_signatures;
        int signature_idx = idx % num_signatures;

        // Get the sample and signature sequences
        int sample_len = d_sample_lengths[sample_idx];
        int signature_len = d_signature_lengths[signature_idx];
        const char *sample_seq = d_sample_sequences[sample_idx];
        const char *sample_qual = d_sample_qualities[sample_idx];
        const char *signature_seq = d_signature_sequences[signature_idx];

        // Initialize the match score (-1 indicates no match)
        double match_score = -1;

        // O(sample_len * signature_len) brute force string matching algorithm
        for (int i = 0; i <= sample_len - signature_len; i++)
        {
            bool match_found = true;
            int match_score_local = 0;
            for (int j = 0; j < signature_len; j++)
            {
                char sample_char = sample_seq[i + j];
                char signature_char = signature_seq[j];
                if (sample_char != signature_char && sample_char != 'N' && signature_char != 'N')
                {
                    match_found = false;
                    break;
                }
                char sample_qual_char = sample_qual[i + j];
                match_score_local += (sample_qual_char - 33);
            }
            if (match_found)
            {
                match_score = (double)match_score_local / signature_len;
                break;
            }
        }

        // Store the match score in the match matrix
        d_match_matrix[sample_idx * num_signatures + signature_idx] = match_score;
    }
}

void runMatcher(
    const std::vector<klibpp::KSeq> &samples, 
    const std::vector<klibpp::KSeq> &signatures, 
    std::vector<MatchResult> &matches)
{
    int num_samples = samples.size();
    int num_signatures = signatures.size();

    // Allocate host arrays for lengths
    int *h_sample_lengths = new int[num_samples];
    int *h_signature_lengths = new int[num_signatures];

    // Allocate char* arrays for sample and signature
    char **h_sample_sequences = new char *[num_samples];
    char **h_signature_sequences = new char *[num_signatures];
    char **h_sample_qualities = new char *[num_samples];

    // Copy data from vectors into arrays
    for (int i = 0; i < num_samples; i++)
    {
        h_sample_lengths[i] = samples[i].seq.size();
        h_sample_sequences[i] = new char[h_sample_lengths[i]];
        h_sample_qualities[i] = new char[h_sample_lengths[i]];
        memcpy(h_sample_sequences[i], samples[i].seq.c_str(), h_sample_lengths[i]);
        memcpy(h_sample_qualities[i], samples[i].qual.c_str(), h_sample_lengths[i]);
    }
    for (int i = 0; i < num_signatures; i++)
    {
        h_signature_lengths[i] = signatures[i].seq.size();
        h_signature_sequences[i] = new char[h_signature_lengths[i]];
        memcpy(h_signature_sequences[i], signatures[i].seq.c_str(), h_signature_lengths[i]);
    }

    // Allocate device memory
    char **d_sample_sequences, **d_sample_qualities;
    char **d_signature_sequences;
    int *d_sample_lengths, *d_signature_lengths;
    double *d_match_matrix;

    cudaMalloc(&d_sample_sequences, num_samples * sizeof(char *));
    cudaMalloc(&d_sample_qualities, num_samples * sizeof(char *));
    cudaMalloc(&d_signature_sequences, num_signatures * sizeof(char *));
    cudaMalloc(&d_sample_lengths, num_samples * sizeof(int));
    cudaMalloc(&d_signature_lengths, num_signatures * sizeof(int));
    cudaMalloc(&d_match_matrix, num_samples * num_signatures * sizeof(double));

    // Copy data to device
    cudaMemcpy(d_sample_lengths, h_sample_lengths, num_samples * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_signature_lengths, h_signature_lengths, num_signatures * sizeof(int), cudaMemcpyHostToDevice);

    // Allocate space and copy data for individual sequences on the GPU
    for (int i = 0; i < num_samples; i++)
    {
        char *d_sample;
        cudaMalloc(&d_sample, h_sample_lengths[i] * sizeof(char));
        cudaMemcpy(d_sample, h_sample_sequences[i], h_sample_lengths[i] * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(&d_sample_sequences[i], &d_sample, sizeof(char *), cudaMemcpyHostToDevice);
        char *d_sample_quality;
        cudaMalloc(&d_sample_quality, h_sample_lengths[i] * sizeof(char));
        cudaMemcpy(d_sample_quality, h_sample_qualities[i], h_sample_lengths[i] * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(&d_sample_qualities[i], &d_sample_quality, sizeof(char *), cudaMemcpyHostToDevice);
    }
    for (int i = 0; i < num_signatures; i++)
    {
        char *d_signature;
        cudaMalloc(&d_signature, h_signature_lengths[i] * sizeof(char));
        cudaMemcpy(d_signature, h_signature_sequences[i], h_signature_lengths[i] * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(&d_signature_sequences[i], &d_signature, sizeof(char *), cudaMemcpyHostToDevice);
    }

    // Configure kernel launch parameters
    int total_pairs = num_samples * num_signatures;
    int threadsPerBlock = 256; // through trial and error
    int blocksPerGrid = (total_pairs + threadsPerBlock - 1) / threadsPerBlock;

    // Launch the string matching kernel
    matchStringsKernel<<<blocksPerGrid, threadsPerBlock>>>(
        d_sample_sequences,
        d_sample_qualities,
        d_sample_lengths,
        d_signature_sequences,
        d_signature_lengths,
        num_samples,
        num_signatures,
        d_match_matrix);
    cudaDeviceSynchronize();

    // Copy match matrix back to host
    double *h_match_matrix = new double[num_samples * num_signatures];
    cudaMemcpy(h_match_matrix, d_match_matrix, num_samples * num_signatures * sizeof(double), cudaMemcpyDeviceToHost);

    // Process the match results
    for (int i = 0; i < num_samples; i++)
    {
        for (int j = 0; j < num_signatures; j++)
        {
            if (h_match_matrix[i * num_signatures + j] != -1)
            {
                MatchResult result = {samples[i].name, signatures[j].name, h_match_matrix[i * num_signatures + j]};
                matches.emplace_back(result);
            }
        }
    }

    // Clean up host and device memory
    delete[] h_sample_lengths;
    delete[] h_signature_lengths;
    for (int i = 0; i < num_samples; i++)
    {
        delete[] h_sample_sequences[i];
        delete[] h_sample_qualities[i];
    }
    delete[] h_sample_sequences;
    delete[] h_sample_qualities;
    for (int i = 0; i < num_signatures; i++)
    {
        delete[] h_signature_sequences[i];
    }
    delete[] h_signature_sequences;
    delete[] h_match_matrix;

    cudaFree(d_sample_sequences);
    cudaFree(d_sample_qualities);
    cudaFree(d_sample_lengths);
    cudaFree(d_signature_sequences);
    cudaFree(d_signature_lengths);
    cudaFree(d_match_matrix);
}
