#include "kseq/kseq.h"
#include "common.h"

// CUDA kernel to perform string matching
__global__ void matchStringsKernel(char **d_samples, int *sample_lengths,
                                   char **d_signatures, int *signature_lengths,
                                   char **d_sample_qualities,
                                   int num_samples, int num_signatures, double *d_match_matrix)
{
    int sample_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int signature_idx = blockIdx.y * blockDim.y + threadIdx.y;

    if (sample_idx < num_samples && signature_idx < num_signatures)
    {
        // Simple string match: Check if the sequences are identical
        int sample_len = sample_lengths[sample_idx];
        int signature_len = signature_lengths[signature_idx];

        // Simple string matching (first occurrence)
        double match_score = -1;
        for (int i = 0; i <= sample_len - signature_len; i++)
        {
            bool match_found = true;
            double match_score_local = 0;
            for (int j = 0; j < signature_len; j++)
            {
                if (d_samples[sample_idx][i + j] != d_signatures[signature_idx][j])
                {
                    match_found = false;
                    break;
                }
                match_score_local += (d_sample_qualities[sample_idx][i + j] - 33);
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

void runMatcher(const std::vector<klibpp::KSeq> &samples, const std::vector<klibpp::KSeq> &signatures, std::vector<MatchResult> &matches)
{
    int num_samples = samples.size();
    int num_signatures = signatures.size();

    // Allocate host arrays for sample and signature lengths
    int *h_sample_lengths = new int[num_samples];
    int *h_signature_lengths = new int[num_signatures];

    // Allocate arrays for sample and signature sequences on the host
    char **h_sample_sequences = new char *[num_samples];
    char **h_signature_sequences = new char *[num_signatures];

    // Allocate arrays for sample qualities on the host
    char **h_sample_qualities = new char *[num_samples];

    // Prepare data to copy to GPU
    for (int i = 0; i < num_samples; ++i)
    {
        h_sample_lengths[i] = samples[i].seq.size();
        h_sample_sequences[i] = new char[h_sample_lengths[i]];
        memcpy(h_sample_sequences[i], samples[i].seq.c_str(), h_sample_lengths[i]);
        h_sample_qualities[i] = new char[h_sample_lengths[i]];
        memcpy(h_sample_qualities[i], samples[i].qual.c_str(), h_sample_lengths[i]);
    }

    for (int i = 0; i < num_signatures; ++i)
    {
        h_signature_lengths[i] = signatures[i].seq.size();
        h_signature_sequences[i] = new char[h_signature_lengths[i]];
        memcpy(h_signature_sequences[i], signatures[i].seq.c_str(), h_signature_lengths[i]);
    }

    // Allocate device memory for sample sequences, signature sequences, and match matrix
    char **d_samples, **d_signatures;
    char **d_sample_qualities;
    int *d_sample_lengths, *d_signature_lengths;
    double *d_match_matrix;

    // Allocate space on GPU for sequences and lengths
    cudaMalloc(&d_samples, num_samples * sizeof(char *));
    cudaMalloc(&d_signatures, num_signatures * sizeof(char *));
    cudaMalloc(&d_sample_qualities, num_samples * sizeof(char *));
    cudaMalloc(&d_sample_lengths, num_samples * sizeof(int));
    cudaMalloc(&d_signature_lengths, num_signatures * sizeof(int));
    cudaMalloc(&d_match_matrix, num_samples * num_signatures * sizeof(double));

    // Copy lengths to GPU
    cudaMemcpy(d_sample_lengths, h_sample_lengths, num_samples * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_signature_lengths, h_signature_lengths, num_signatures * sizeof(int), cudaMemcpyHostToDevice);

    // Allocate space for individual sequences on the GPU
    for (int i = 0; i < num_samples; ++i)
    {
        char *d_sample;
        cudaMalloc(&d_sample, h_sample_lengths[i] * sizeof(char));
        cudaMemcpy(d_sample, h_sample_sequences[i], h_sample_lengths[i] * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(&d_samples[i], &d_sample, sizeof(char *), cudaMemcpyHostToDevice);

        char *d_sample_quality;
        cudaMalloc(&d_sample_quality, h_sample_lengths[i] * sizeof(char));
        cudaMemcpy(d_sample_quality, h_sample_qualities[i], h_sample_lengths[i] * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(&d_sample_qualities[i], &d_sample_quality, sizeof(char *), cudaMemcpyHostToDevice);
    }

    for (int i = 0; i < num_signatures; ++i)
    {
        char *d_signature;
        cudaMalloc(&d_signature, h_signature_lengths[i] * sizeof(char));
        cudaMemcpy(d_signature, h_signature_sequences[i], h_signature_lengths[i] * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(&d_signatures[i], &d_signature, sizeof(char *), cudaMemcpyHostToDevice);
    }

    // Configure kernel launch parameters
    dim3 threadsPerBlock(16, 16);
    dim3 blocksPerGrid((num_samples + threadsPerBlock.x - 1) / threadsPerBlock.x,
                       (num_signatures + threadsPerBlock.y - 1) / threadsPerBlock.y);

    // Launch the string matching kernel
    matchStringsKernel<<<blocksPerGrid, threadsPerBlock>>>(d_samples, d_sample_lengths,
                                                           d_signatures, d_signature_lengths,
                                                           d_sample_qualities,
                                                           num_samples, num_signatures, d_match_matrix);

    // Copy match matrix back to host
    double *h_match_matrix = new double[num_samples * num_signatures];
    cudaMemcpy(h_match_matrix, d_match_matrix, num_samples * num_signatures * sizeof(double), cudaMemcpyDeviceToHost);

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
    delete[] h_match_matrix;

    for (int i = 0; i < num_samples; ++i)
    {
        delete[] h_sample_sequences[i];
    }

    for (int i = 0; i < num_signatures; ++i)
    {
        delete[] h_signature_sequences[i];
    }

    delete[] h_sample_sequences;
    delete[] h_signature_sequences;

    cudaFree(d_samples);
    cudaFree(d_signatures);
    cudaFree(d_sample_lengths);
    cudaFree(d_signature_lengths);
    cudaFree(d_match_matrix);
}
