#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>

std::vector<std::string> compute_prod (int N) {
    if (N <= 0)
        return std::vector<std::string> (0);
    std::vector<std::string> prod ((1 << N));
    for (int ii = 0; ii < (1 << N); ++ii) {
        prod[ii].reserve(N);
        for (int jj = 0; jj < N; ++jj)
            prod[ii] += '0' + ((ii >> (N - jj - 1)) & 1);
    }
    return prod;
}

void precompute_products (int N, std::unordered_map<int, std::vector<std::string>> &prods) {
    for (int ii = -1; ii <= (N + 2); ++ii)
        prods[ii] = std::move(compute_prod(ii));
}

std::vector<std::vector<int>> transpose(const std::vector<std::vector<int>> &X, int start, int end) {
    // Compute the dimensions of the transposed vector
    int tlen = X[start].size();
    if (end == -1)
        end = X.size();
    std::vector<std::vector<int>> t(tlen, std::vector<int>(end - start));

    // Fill in the transposed vector
    for (int ii = start; ii < end; ++ii)
        for (int jj = 0; jj < tlen; ++jj)
            t[jj][ii - start] = X[ii][jj];
    return t;
}

// Write the phased haplotypes to a file
void print_haplotypes(const std::vector<std::vector<int>> &h_t, std::string outf) {
    std::ofstream of (outf, std::ios::app);
    of.seekp(0, std::ios::end);
    for (std::vector<int> row : h_t) {
        int rlen = row.size();
        for (int ii = 0; ii < rlen; ++ii) {
            if (ii)
                of << " ";
            of << row[ii];
        }
        of << "\n";
    }
    of.close();
}


std::string complement(const std::string &perm) {
    int hlen = perm.size();
    std::string c;
    c.reserve(hlen);
    for (int ii = 0; ii < hlen; ++ii)
        c += perm[ii] == '1' ? '0' : '1';
    return c;
}

std::string gen_hapl(const std::vector<int> &geno, std::string perm, std::string mperm) {
    int pidx = 0, midx = 0;
    int glen = geno.size();
    std::string hapl;
    hapl.reserve(glen);
    for (int ii = 0; ii < glen; ++ii) {
        if (geno[ii] == 0)
            hapl += '0';
        else if (geno[ii] == 2)
            hapl += '1';
        else if (geno[ii] == -1) {
            hapl += mperm[midx];
            midx++;
        }
        else {
            hapl += perm[pidx];
            pidx++;
        }
    }
    return hapl;
}

void generate_and_initialize(   const std::vector<std::vector<int>> & genos,
                                std::vector<std::vector<std::tuple<std::string,std::string,float>>> &compat_probs,
                                std::unordered_map<std::string, float> &hapl_probs,
                                std::unordered_map<int, std::vector<std::string>> &prods)
{
    int glen = genos.size();
    for (int ii = 0; ii < glen; ii++) {
        int hetero_cnt = std::count(genos[ii].begin(), genos[ii].end(), 1);
        int missing_cnt = std::count(genos[ii].begin(), genos[ii].end(), -1);
        int hposs = hetero_cnt ? 1 << (hetero_cnt - 1) : 1;
        int mposs = 1 << missing_cnt;
        int total_poss = hposs * mposs;
        compat_probs[ii] = std::vector<std::tuple<std::string,std::string,float>> (total_poss);
        if (!hetero_cnt) {
            if (!missing_cnt) {
                std::string homo_hapl = gen_hapl(genos[ii], "", "");
                compat_probs[ii][0] = std::make_tuple(homo_hapl, homo_hapl, 1.0 / total_poss);
                hapl_probs[homo_hapl] = 0.0;
                continue;
            }
            else {
                for (int kk = 0; kk < (1 << missing_cnt); ++kk) {
                    std::string homo_hapl = gen_hapl(genos[ii], "", prods[missing_cnt][kk]);
                    compat_probs[ii][kk] = std::make_tuple(homo_hapl, homo_hapl, 1.0 / total_poss);
                    hapl_probs[homo_hapl] = 0.0;
                }
                continue;
            }
        }
        for (int kk = 0; kk < hposs; ++kk) {
            if (!missing_cnt) {
                std::string hapl = gen_hapl(genos[ii], prods[hetero_cnt][kk], "");
                std::string c_hapl = gen_hapl(genos[ii], complement(prods[hetero_cnt][kk]), "");
                compat_probs[ii][kk] = std::make_tuple(hapl, c_hapl, 1.0 / total_poss);
                hapl_probs[hapl] = 0.0;
                hapl_probs[c_hapl] = 0.0;
            }
            else {
                for (int ll = 0; ll < mposs; ++ll) {
                    std::string hapl = gen_hapl(genos[ii], prods[hetero_cnt][kk], prods[missing_cnt][ll]);
                    std::string c_hapl = gen_hapl(genos[ii], complement(prods[hetero_cnt][kk]), prods[missing_cnt][ll]);
                    compat_probs[ii][kk * mposs + ll] = std::make_tuple(hapl, c_hapl, 1.0 / total_poss);
                    hapl_probs[hapl] = 0.0;
                    hapl_probs[c_hapl] = 0.0;
                }
                continue;
            }
        }
    }
}

// This version of Hamming distance is dangerous. Namely, h1 might well be
// longer than h0.
int hd (const std::vector<int> &h0, const std::vector<int> &h1) {
    int hlen = h0.size();
    int dist = 0;
    for (int ii = 0; ii < hlen; ii++)
        if (h0[ii] != h1[ii])
            dist++;
    return dist;
}

void check_for_parity (std::vector<std::vector<int>> &hapl, const std::vector<std::vector<int>> end_bits) {
    int num_hapls = hapl.size();
    for (int ii = 0; ii < num_hapls; ii += 2) {
        // Yes. I know hapl[ii] and hapl[ii + 1] are longer than end_bits.
        if (hd(end_bits[ii], hapl[ii + 1]) < hd(end_bits[ii], hapl[ii]))
            std::swap(hapl[ii], hapl[ii + 1]);
    }
}

void me_step(   std::vector<std::vector<std::tuple<std::string,std::string,float>>> &compat_probs,
                std::unordered_map<std::string, float> &hapl_probs)
{
    // Set all values in hapl_probs to 0.0
    for (std::unordered_map<std::string, float>::iterator hp = hapl_probs.begin(); hp != hapl_probs.end(); ++hp)
        hapl_probs[hp->first] = 0.0;

    // M step
    int clen = compat_probs.size();
    #pragma omp parallel for
    for (int ii = 0; ii < clen; ii++) {
        int ilen = compat_probs[ii].size();
        for (int jj = 0; jj < ilen; jj++) {
            float cprob = std::get<2>(compat_probs[ii][jj]);
            hapl_probs[std::get<0>(compat_probs[ii][jj])] += cprob / (2.0 * clen);
            hapl_probs[std::get<1>(compat_probs[ii][jj])] += cprob / (2.0 * clen);
        }
    }

    // E step
    #pragma omp parallel for
    for (int ii = 0; ii < clen; ii++) {
        int ilen = compat_probs[ii].size();
        float ph_sum = 0.0;
        for (int jj = 0; jj < ilen; jj++) {
            float ph1_ph2 = hapl_probs[std::get<0>(compat_probs[ii][jj])] * hapl_probs[std::get<1>(compat_probs[ii][jj])];
            std::get<2>(compat_probs[ii][jj]) = ph1_ph2;
            ph_sum += ph1_ph2;
        }
        for (int jj = 0; jj < ilen; jj++) {
            std::get<2>(compat_probs[ii][jj]) /= ph_sum;
        }
    }

    return;
}

float q_metric( std::vector<std::vector<std::tuple<std::string,std::string,float>>> &compat_probs,
                std::unordered_map<std::string, float> &hapl_probs)
{
    float q_sum = 0.0;
    int clen = compat_probs.size();
    for (int ii = 0; ii < clen; ii++) {
        float max_p = std::get<2>(compat_probs[ii][0]);
        int max_idx = 0;
        int num_compat = compat_probs[ii].size();
        for (int jj = 1; jj < num_compat; jj++) {
            float p = std::get<2>(compat_probs[ii][jj]);
            if (p > max_p) {
                max_p = p;
                max_idx = jj;
            }
        }

        q_sum += log(hapl_probs[std::get<0>(compat_probs[ii][max_idx])]);
        q_sum += log(hapl_probs[std::get<1>(compat_probs[ii][max_idx])]);
    }
    return q_sum;
}

std::vector<int> to_vector(std::string h) {
    int hlen = h.size();
    std::vector<int> v(hlen);
    for (int ii = 0; ii < hlen; ii++)
        v[ii] = h[ii] - '0';
    return v;
}

void process_EM(std::vector<std::vector<int>> genos, std::string of, std::vector<std::vector<int>> &end_bits, int overlap, std::unordered_map<int, std::vector<std::string>> &prods) {
    // For each individual, this contains pairs of compatible haplotypes
    std::vector<std::vector<std::tuple<std::string,std::string,float>>> compat_probs(genos.size(), std::vector<std::tuple<std::string, std::string, float>>(0));
    // Not 100% sure this is sufficient precision. May need to switch to double.
    std::unordered_map<std::string,float> hapl_probs;
    generate_and_initialize(genos, compat_probs, hapl_probs, prods);
    int glen = genos.size();
    int window = genos[0].size();
    for (int ii = 0; ii < 6; ii++) {
        me_step(compat_probs, hapl_probs);
        //std::cout << q_metric(compat_probs, hapl_probs) << std::endl;
    }

    std::vector<std::vector<int>> inferred_hapls(2 * glen, std::vector<int>(window));
    for (int ii = 0; ii < glen; ii++) {
        float max_p = std::get<2>(compat_probs[ii][0]);
        int max_idx = 0;
        int num_compat = compat_probs[ii].size();
        for (int jj = 1; jj < num_compat; jj++) {
            float p = std::get<2>(compat_probs[ii][jj]);
            if (p > max_p) {
                max_p = p;
                max_idx = jj;
            }
        }

        inferred_hapls[2 * ii] = to_vector(std::get<0>(compat_probs[ii][max_idx]));
        inferred_hapls[2 * ii + 1] = to_vector(std::get<1>(compat_probs[ii][max_idx]));
    }

    if (end_bits.size())
        check_for_parity(inferred_hapls, end_bits);
    std::vector<std::vector<int>> h_t = transpose(inferred_hapls, 0, -1);
    if (end_bits.size()) {
        std::vector<std::vector<int>> h_slice (h_t.begin() + overlap, h_t.end());
        print_haplotypes(h_slice, of);
    }
    else
        print_haplotypes(h_t, of);
    std::vector<std::vector<int>> end_slice (h_t.end() - overlap, h_t.end());
    end_bits = transpose(end_slice, 0, -1);
}

void general_stuff(int t_no, int t_cnt, int num_snps, int num_indivs, int hm_max, int overlap, std::string pfx) {
    std::unordered_map<int, std::vector<std::string>> prods;
    precompute_products(hm_max, prods);

    // Need to figure out the start and end of this portion
    int start = t_no * (int) (num_snps / t_cnt);
    int maybe_end =  (t_no + 1) * (int) (num_snps / t_cnt) + overlap;
    int end = t_no != (t_cnt - 1) ? maybe_end : num_snps;

    // Make and fill a vector of the requisite dimensions
    std::vector<std::vector<int>> X(end - start, std::vector<int>(num_indivs));
    int row_idx = 0;
    std::string masked_f = pfx + "_masked.txt";
    std::string of = pfx + "_sol_wip0_" + std::to_string(t_no) + ".txt";
    truncate(of.c_str(), 0);

    std::string line;
    std::ifstream mf0(masked_f);
    while(std::getline(mf0, line)) {
        if (row_idx >= start && row_idx < end) {
            std::istringstream is(line);
            std::string g;
            int col_idx = 0;
            while(getline(is, g, ' ')) {
                g = !g.compare("*") ? "-1" : g;
                X[row_idx - start][col_idx] = std::stoi(g);
                col_idx++;
            }
        }
        row_idx++;
    }

    std::vector<std::vector<int>> end_bits;
    int idx = 0;
    std::vector<int> num_hm(num_indivs);
    for (int ii = 0; ii < (end - start); ++ii) {
        for (int jj = 0; jj < num_indivs; ++jj)
            if (X[ii][jj] == 1 || X[ii][jj] == -1)
                num_hm[jj]++;

        if (*max_element(num_hm.begin(), num_hm.end()) >= hm_max) {
            // No point in leaving the last couple SNPs lonely
            if (ii < end - start - 2) {
                process_EM(std::move(transpose(X, idx, ii)), of, end_bits, overlap, prods);
                idx = ii - overlap;
                std::fill(num_hm.begin(), num_hm.end(), 0);
                // Need to fill num_hm within the overlapped range
                for (int jj = idx; jj <= ii; ++jj)
                    for (int kk = 0; kk < num_indivs; ++kk)
                        if (X[jj][kk] == 1 || X[jj][kk] == -1)
                            ++num_hm[kk];
                continue;
            }
        }
        if (ii == end - start - 1)
            process_EM(std::move(transpose(X, idx, -1)), of, end_bits, overlap, prods);
    }
}

int main(int argc, char* argv[])
{
    // Very basic argument checking
    if (argc != 5) {
        std::cerr << "Please provide a file prefix, thread count, heterozygous count, and overlap" << std::endl;
        return 1;
    }

    const int thread_cnt = std::stoi(argv[2]);
    // Sum of heterozygous and missing. Max for an individual in a window.
    const int hm_max = std::stoi(argv[3]);
    // Overlap hyperparam: haplotype with smaller Hamming distance will be chosen
    const int overlap = std::stoi(argv[4]);

    // Generate filenames and initialize file parsing-related variables
    std::string masked_f = std::string(argv[1]) + "_masked.txt";
    std::string output_f = std::string(argv[1]) + "_sol_wip0.txt";
    std::cout << masked_f << std::endl << output_f << std::endl;
    truncate(output_f.c_str(), 0);
    int num_snps = 0, num_indivs = 0;
    std::string line;
    std::ifstream mf(masked_f);

    // Count number of SNPs and individuals
    while(std::getline(mf, line)) {
        if (!num_snps) {
            std::istringstream is(line);
            std::string g;
            while(getline(is, g, ' '))
                num_indivs++;
        }
        num_snps++;
    }
    mf.close();
    std::cout << "SNPs: " << num_snps << std::endl;
    std::cout << "Individuals: " << num_indivs << std::endl;

    std::vector<std::thread> threads;
    for (int ii = 0; ii < thread_cnt; ii++)
        threads.push_back(std::thread (general_stuff, ii, thread_cnt, num_snps, num_indivs, hm_max, overlap, std::string(argv[1])));

    std::vector<std::vector<int>> end_bits;
    // Combine the threads' results
    for (int ii = 0; ii < thread_cnt; ii++) {
        threads[ii].join();
        // Need to figure out the start and end of this portion
        int start = ii * (int) (num_snps / thread_cnt);
        int maybe_end =  (ii + 1) * (int) (num_snps / thread_cnt) + overlap;
        int end = ii != (thread_cnt - 1) ? maybe_end : num_snps;

        std::cout << ii << ": " << start << "," << end << std::endl;

        std::vector<std::vector<int>> Xii(end - start, std::vector<int>(2 * num_indivs));
        // Read the thread's output file
        std::string of = std::string(argv[1]) + "_sol_wip0_" + std::to_string(ii) + ".txt";
        std::string line;
        std::ifstream xf(of);
        int row_idx = 0;

        // Count number of SNPs and individuals
        while(std::getline(xf, line)) {
            std::istringstream is(line);
            std::string g;
            int col_idx = 0;
            while(getline(is, g, ' ')) {
                Xii[row_idx][col_idx] = std::stoi(g);
                col_idx++;
            }
            row_idx++;
        }

        // This part is terribly wasteful. Feel free to make a version
        // of check_for_parity that works on the transposed matrix.
        std::vector<std::vector<int>> inferred_hapls = transpose(Xii, 0, -1);
        if (end_bits.size())
            check_for_parity(inferred_hapls, end_bits);
        std::vector<std::vector<int>> h_t = transpose(inferred_hapls, 0, -1);
        if (end_bits.size()) {
            std::vector<std::vector<int>> h_slice (h_t.begin() + overlap, h_t.end());
            print_haplotypes(h_slice, output_f);
        }
        else
            print_haplotypes(h_t, output_f);
        std::vector<std::vector<int>> end_slice (Xii.end() - overlap, Xii.end());
        end_bits = transpose(end_slice, 0, -1);
    }
}
