#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <omp.h>

using Eigen::Map;
using Eigen::MatrixXd;
using namespace Rcpp;

struct Sorter
{
    int index;
    double value;
};

bool operator<(const Sorter &lhs, const Sorter &rhs)
{
    return lhs.value > rhs.value; //after sort vector will be big to small!
}

// [[Rcpp::export]]
SEXP DsGene(SEXP expr_r, SEXP label_r, bool as_df, int num_top, int additional, int n_threads)
{
    if (n_threads)
        omp_set_num_threads(n_threads);
    int n_cell = LENGTH(VECTOR_ELT(expr_r, 0));
    int n_gene = LENGTH(expr_r) - 1;
    std::vector<int> label = as<std::vector<int> >(label_r); // label_r is 0:n_cell - 1
    int n_label = *std::max_element(label.begin(), label.end()) + 1;
    std::vector<int> label_n(n_label);
    for (int i = 0; i < n_cell; ++i)
        label_n[label[i]]++;
    std::vector<std::vector<double> *> type_expr; // row is gene, col is label
    for (int i = 0; i < n_gene; ++i)
        type_expr.push_back(new std::vector<double>(n_label));
    std::vector<Sorter> sum_ds(n_gene);
    #pragma omp parallel for simd
    for (int i = 0; i < n_gene; ++i)
    {
        std::vector<double> *cur_row = type_expr[i];
        for (int j = 0; j < n_cell; ++j)
            cur_row->at(label[j]) += REAL(VECTOR_ELT(expr_r, i))[j];
        double sum = 0;
        double s_ds = 0;
        for (int j = 0; j < n_label; ++j)
        {
            cur_row->at(j) /= label_n[j];
            sum += cur_row->at(j);
            cur_row->at(j) = std::log(cur_row->at(j) + 1);
            s_ds += cur_row->at(j);
        }
        sum_ds[i] = {i, std::log(sum / n_label + 1) * n_label - s_ds};
    }
    if (as_df)
    {
        std::vector<double> ret_data;
        ret_data.reserve(n_gene * (n_label + 1));
        for (int i = 0; i < n_gene; ++i)
        {
            std::copy(type_expr[i]->begin(), type_expr[i]->end(), std::back_inserter(ret_data));
            ret_data.push_back(sum_ds[i].value);
        }
        #pragma omp simd
        for (int i = 0; i < n_gene; ++i)
        {
            delete(type_expr[i]);
            type_expr[i] = NULL;
        }
        NumericMatrix ret(n_label + 1, n_gene, ret_data.begin());
        return transpose(ret);
    }
    std::vector<int> ind(num_top);
    if (num_top)
    {
        std::sort(sum_ds.begin(), sum_ds.end());
        int i = 0;
        for (int k = 0; k < num_top; ++i, ++k)
            ind[k] = sum_ds[i].index;
    }
    if (additional)
    {
        std::vector<int> addi(n_label * additional);
        #pragma omp parallel for
        for (int j = 0; j < n_label; ++j)
        {
            std::vector<Sorter> sorted(n_gene);
            for (int i = 0; i < n_gene; ++i)
                sorted[i] = {i, type_expr[i]->at(j)};
            std::sort(sorted.begin(), sorted.end());
            int i = 0;
            for (int k = 0; k < additional; ++i, ++k)
                addi[j * additional + k] = sorted[i].index;
        }
        std::copy(addi.begin(), addi.end(), std::back_inserter(ind));
        std::sort(ind.begin(), ind.end());
        ind.erase(std::unique(ind.begin(), ind.end()), ind.end());
    }
    #pragma omp simd
    for (int i = 0; i < n_gene; ++i)
    {
        delete(type_expr[i]);
        type_expr[i] = NULL;
    }
    IntegerVector ind_out = wrap(ind);
    return ind_out;
}

// [[Rcpp::export]]
SEXP GenProb(SEXP expr_r, SEXP label_r, SEXP geneset_r, int n_threads)
{
    if (n_threads)
        omp_set_num_threads(n_threads);
    std::vector<int> label = as<std::vector<int> >(label_r); // label_r is 0:n_cell - 1
    std::vector<int> genes = as<std::vector<int> >(geneset_r);
    int n_cell = LENGTH(VECTOR_ELT(expr_r, 0));
    int n_gene = genes.size();
    int n_label = *std::max_element(label.begin(), label.end()) + 1;
    std::vector<double> prob(n_label * n_gene);
    #pragma omp parallel for simd
    for (int i = 0; i < n_gene; ++i)
        for (int j = 0; j < n_cell; ++j)
            prob[label[j] * n_gene + i] += std::log(REAL(VECTOR_ELT(expr_r, genes[i]))[j] + 1) / std::log(2.0);
    std::vector<double> sums(n_label);
    for (int i = 0; i < n_gene; ++i)
        for (int j = 0; j < n_label; ++j)
            sums[j] += prob[j * n_gene + i];
    #pragma omp parallel for simd
    for (int i = 0; i < n_gene; ++i)
        for (int j = 0; j < n_label; ++j)
            prob[j * n_gene + i] = std::log(prob[j * n_gene + i] + 1) - std::log(sums[j] + n_gene);
    NumericMatrix ret(n_gene, n_label, prob.begin());
    return ret;
}

// [[Rcpp::export]]
SEXP Gambler(SEXP expr_r, SEXP prob_r, bool ret_tab, int n_threads)
{
    if (n_threads)
        omp_set_num_threads(n_threads);
    const Map<MatrixXd> expr(as<Map<MatrixXd> >(expr_r));
    const Map<MatrixXd> prob(as<Map<MatrixXd> >(prob_r));
    MatrixXd ret = expr * prob;
    int m = ret.rows();
    int n = ret.cols();
    std::vector<int> ind(m);
    #pragma omp parallel for simd
    for (int i = 0; i < m; ++i)
    {
        double v = ret(i, 0);
        for (int j = 0; j < n; ++j)
        {
            double u = ret(i, j);
            if (u > v)
            {
                v = u;
                ind[i] = j;
            }
        }
        if (ret_tab)
        {
            double sum = 0;
            for (int j = 0; j < n; ++j)
            {
                double a = std::exp(ret(i, j) - v);
                ret(i, j) = a;
                sum += a;
            }
            for (int j = 0; j < n; ++j)
                ret(i, j) /= sum;
        }
    }
    if (ret_tab)
      return wrap(ret);
    IntegerVector ind_out = wrap(ind);
    return ind_out;
}

// [[Rcpp::export]]
SEXP GenEntr(SEXP expr_r, double window, int n_threads)
{
    if (n_threads)
        omp_set_num_threads(n_threads);
    int n_cell = LENGTH(VECTOR_ELT(expr_r, 0));
    int n_gene = LENGTH(expr_r);
    std::vector<double> entropy(n_gene);
    #pragma omp parallel for simd
    for (int i = 0; i < n_gene; ++i)
    {
        double sum = 0;
        int states = std::ceil(1000000.0 / window);
        std::vector<double> discretize(states);
        for (int j = 0; j < n_cell; ++j)
            discretize[std::ceil(REAL(VECTOR_ELT(expr_r, i))[j] / window)]++;
        for (int j = 0; j < states; ++j)
            sum += discretize[j];
        for (int j = 0; j < states; ++j)
            if (discretize[j])
                entropy[i] -= discretize[j] / sum * std::log(discretize[j] / sum);
    }
    NumericVector entr = wrap(entropy);
    return entr;
}

// [[Rcpp::export]]
SEXP NullTest(SEXP ref_r, SEXP query_r, SEXP null_r, SEXP label_r, int n_feat, int n_threads)
{
    if (n_threads)
        omp_set_num_threads(n_threads);
    const Map<MatrixXd> ref(as<Map<MatrixXd> >(ref_r));
    const Map<MatrixXd> query(as<Map<MatrixXd> >(query_r));
    int n_ref = ref.rows();
    int n_query = query.rows();
    int n_gene = ref.cols();
    std::vector<int> label = as<std::vector<int> >(label_r); // label_r is 0:n_cell - 1
    std::vector<double> null = as<std::vector<double> >(null_r);
    int n_label = *std::max_element(label.begin(), label.end()) + 1;
    std::vector<int> label_n(n_label);
    for (int i = 0; i < n_ref; ++i)
        label_n[label[i]]++;
    std::vector<std::vector<double> *> type_expr; // row is gene, col is label
    for (int i = 0; i < n_gene; ++i)
        type_expr.push_back(new std::vector<double>(n_label));
    std::vector<double> e1(n_gene);
    std::vector<Sorter> ds(n_gene);
    #pragma omp parallel for simd
    for (int i = 0; i < n_gene; ++i)
    {
        std::vector<double> *cur_row = type_expr[i];
        for (int j = 0; j < n_ref; ++j)
            cur_row->at(label[j]) += ref(j, i);
        double sum = 0;
        for (int j = 0; j < n_label; ++j)
        {
            cur_row->at(j) /= label_n[j];
            sum += cur_row->at(j);
        }
        e1[i] = sum / n_label;
        ds[i] = {i, std::log(e1[i] + 1) - std::log(null[i] + 1)};
    }
    double sum_e1 = 0;
    double sum_null = 0;
    for (int i = 0; i < n_gene; ++i)
    {
        sum_e1 += e1[i];
        sum_null += null[i];
    }
    #pragma omp simd
    for (int i = 0; i < n_gene; ++i)
        e1[i] = std::log((e1[i] + 1) / (sum_e1 + n_gene)) - std::log((null[i] + 1) / (sum_null + n_gene));
    std::sort(ds.begin(), ds.end());
    std::vector<double> prob(n_query);
    #pragma omp parallel for simd
    for (int j = 0; j < n_query; ++j)
        for (int i = 0; i < n_feat; ++i)
            prob[j] += e1[ds[i].index] * query(j, ds[i].index);
    #pragma omp simd
    for (int i = 0; i < n_gene; ++i)
    {
        delete(type_expr[i]);
        type_expr[i] = NULL;
    }
    return wrap(prob);
}
