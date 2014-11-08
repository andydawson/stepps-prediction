// manual gradient for the prediction model

#define EIGEN_DONT_PARALLELIZE

#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <typeinfo>
#include <cmath>

#include <boost/exception/all.hpp>

#include <stan/model/prob_grad.hpp>
#include <stan/prob/distributions.hpp>
#include <stan/math/matrix.hpp>
#include <stan/math.hpp>
#include <stan/io/dump.hpp>
#include <stan/io/reader.hpp>
#include <stan/io/writer.hpp>
#include <stan/io/csv_writer.hpp>

#include <unsupported/Eigen/KroneckerProduct>

#include "timer.hpp"

namespace pred_model_namespace {

  using namespace std;
  using namespace stan::math;
  using namespace stan::prob;
  using namespace stan::io;
  using namespace Eigen;

  using stan::math::lgamma;

  typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
  typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
  typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

  class pred_model : public stan::model::prob_grad {
private:
    int K;
    int N;
    int T;
    int N_knots;
    int N_cores;
    vector<vector<int> > y;
    vector_d rho;
    vector_d eta;
    double gamma;
    double psi;
    vector_d phi;
    vector<int> idx_cores;
    matrix_d d;
    matrix_d d_knots;
    matrix_d d_inter;
    matrix_d w;
    matrix_d lag;
    int N_p;
    double P;
    int W;
    vector_d eta2;
    vector_d zeros;
    matrix_d Eye_knots;
    vector<matrix_d> C_s;
    vector<Eigen::LLT<matrix_d> > llt_s;
    vector<matrix_d> C_s_L;
    vector<matrix_d> C_s_inv;
    vector<matrix_d> c_s;
    matrix_d M;

    Eigen::HouseholderQR<vector_d> qr;
    matrix_d fatQ;
    vector_d Q;
    row_vector_d QT;

    static double const ksi_lb = 0;
    static double const ksi_ub = 20;
    static double const sigma_lb = 0;
    static double const sigma_ub = 100;
    static double const lambda_lb = 0;
    static double const lambda_ub = 2;

public:
    pred_model(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad::prob_grad(0) {
        static const char* function__ = "pred_model_namespace::pred_model(%1%)";
        (void) function__; // dummy call to supress warning
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        context__.validate_dims("data initialization", "K", "int", context__.to_vec());
        K = int(0);
        vals_i__ = context__.vals_i("K");
        pos__ = 0;
        K = vals_i__[pos__++];
        context__.validate_dims("data initialization", "N", "int", context__.to_vec());
        N = int(0);
        vals_i__ = context__.vals_i("N");
        pos__ = 0;
        N = vals_i__[pos__++];
        context__.validate_dims("data initialization", "T", "int", context__.to_vec());
        T = int(0);
        vals_i__ = context__.vals_i("T");
        pos__ = 0;
        T = vals_i__[pos__++];
        context__.validate_dims("data initialization", "N_knots", "int", context__.to_vec());
        N_knots = int(0);
        vals_i__ = context__.vals_i("N_knots");
        pos__ = 0;
        N_knots = vals_i__[pos__++];
        context__.validate_dims("data initialization", "N_cores", "int", context__.to_vec());
        N_cores = int(0);
        vals_i__ = context__.vals_i("N_cores");
        pos__ = 0;
        N_cores = vals_i__[pos__++];
        context__.validate_dims("data initialization", "y", "int", context__.to_vec((N_cores * T),K));
        stan::math::validate_non_negative_index("y", "(N_cores * T)", (N_cores * T));
        stan::math::validate_non_negative_index("y", "K", K);
        y = std::vector<std::vector<int> >((N_cores * T),std::vector<int>(K,int(0)));
        vals_i__ = context__.vals_i("y");
        pos__ = 0;
        size_t y_limit_1__ = K;
        for (size_t i_1__ = 0; i_1__ < y_limit_1__; ++i_1__) {
            size_t y_limit_0__ = (N_cores * T);
            for (size_t i_0__ = 0; i_0__ < y_limit_0__; ++i_0__) {
                y[i_0__][i_1__] = vals_i__[pos__++];
            }
        }
        stan::math::validate_non_negative_index("rho", "(K - 1)", (K - 1));
        rho = vector_d((K - 1));
        context__.validate_dims("data initialization", "rho", "vector_d", context__.to_vec((K - 1)));
        vals_r__ = context__.vals_r("rho");
        pos__ = 0;
        size_t rho_i_vec_lim__ = (K - 1);
        for (size_t i_vec__ = 0; i_vec__ < rho_i_vec_lim__; ++i_vec__) {
            rho[i_vec__] = vals_r__[pos__++];
        }
        stan::math::validate_non_negative_index("eta", "(K - 1)", (K - 1));
        eta = vector_d((K - 1));
        context__.validate_dims("data initialization", "eta", "vector_d", context__.to_vec((K - 1)));
        vals_r__ = context__.vals_r("eta");
        pos__ = 0;
        size_t eta_i_vec_lim__ = (K - 1);
        for (size_t i_vec__ = 0; i_vec__ < eta_i_vec_lim__; ++i_vec__) {
            eta[i_vec__] = vals_r__[pos__++];
        }
        context__.validate_dims("data initialization", "gamma", "double", context__.to_vec());
        gamma = double(0);
        vals_r__ = context__.vals_r("gamma");
        pos__ = 0;
        gamma = vals_r__[pos__++];
        context__.validate_dims("data initialization", "psi", "double", context__.to_vec());
        psi = double(0);
        vals_r__ = context__.vals_r("psi");
        pos__ = 0;
        psi = vals_r__[pos__++];
        stan::math::validate_non_negative_index("phi", "K", K);
        phi = vector_d(K);
        context__.validate_dims("data initialization", "phi", "vector_d", context__.to_vec(K));
        vals_r__ = context__.vals_r("phi");
        pos__ = 0;
        size_t phi_i_vec_lim__ = K;
        for (size_t i_vec__ = 0; i_vec__ < phi_i_vec_lim__; ++i_vec__) {
            phi[i_vec__] = vals_r__[pos__++];
        }
        context__.validate_dims("data initialization", "idx_cores", "int", context__.to_vec(N_cores));
        stan::math::validate_non_negative_index("idx_cores", "N_cores", N_cores);
        idx_cores = std::vector<int>(N_cores,int(0));
        vals_i__ = context__.vals_i("idx_cores");
        pos__ = 0;
        size_t idx_cores_limit_0__ = N_cores;
        for (size_t i_0__ = 0; i_0__ < idx_cores_limit_0__; ++i_0__) {
            idx_cores[i_0__] = vals_i__[pos__++];
        }
        context__.validate_dims("data initialization", "d", "matrix_d", context__.to_vec(N,N));
        stan::math::validate_non_negative_index("d", "N", N);
        stan::math::validate_non_negative_index("d", "N", N);
        d = matrix_d(N,N);
        vals_r__ = context__.vals_r("d");
        pos__ = 0;
        size_t d_m_mat_lim__ = N;
        size_t d_n_mat_lim__ = N;
        for (size_t n_mat__ = 0; n_mat__ < d_n_mat_lim__; ++n_mat__) {
            for (size_t m_mat__ = 0; m_mat__ < d_m_mat_lim__; ++m_mat__) {
                d(m_mat__,n_mat__) = vals_r__[pos__++];
            }
        }
        context__.validate_dims("data initialization", "d_knots", "matrix_d", context__.to_vec(N_knots,N_knots));
        stan::math::validate_non_negative_index("d_knots", "N_knots", N_knots);
        stan::math::validate_non_negative_index("d_knots", "N_knots", N_knots);
        d_knots = matrix_d(N_knots,N_knots);
        vals_r__ = context__.vals_r("d_knots");
        pos__ = 0;
        size_t d_knots_m_mat_lim__ = N_knots;
        size_t d_knots_n_mat_lim__ = N_knots;
        for (size_t n_mat__ = 0; n_mat__ < d_knots_n_mat_lim__; ++n_mat__) {
            for (size_t m_mat__ = 0; m_mat__ < d_knots_m_mat_lim__; ++m_mat__) {
                d_knots(m_mat__,n_mat__) = vals_r__[pos__++];
            }
        }
        context__.validate_dims("data initialization", "d_inter", "matrix_d", context__.to_vec(N,N_knots));
        stan::math::validate_non_negative_index("d_inter", "N", N);
        stan::math::validate_non_negative_index("d_inter", "N_knots", N_knots);
        d_inter = matrix_d(N,N_knots);
        vals_r__ = context__.vals_r("d_inter");
        pos__ = 0;
        size_t d_inter_m_mat_lim__ = N;
        size_t d_inter_n_mat_lim__ = N_knots;
        for (size_t n_mat__ = 0; n_mat__ < d_inter_n_mat_lim__; ++n_mat__) {
            for (size_t m_mat__ = 0; m_mat__ < d_inter_m_mat_lim__; ++m_mat__) {
                d_inter(m_mat__,n_mat__) = vals_r__[pos__++];
            }
        }
        context__.validate_dims("data initialization", "w", "matrix_d", context__.to_vec(N_cores,N));
        stan::math::validate_non_negative_index("w", "N_cores", N_cores);
        stan::math::validate_non_negative_index("w", "N", N);
        w = matrix_d(N_cores,N);
        vals_r__ = context__.vals_r("w");
        pos__ = 0;
        size_t w_m_mat_lim__ = N_cores;
        size_t w_n_mat_lim__ = N;
        for (size_t n_mat__ = 0; n_mat__ < w_n_mat_lim__; ++n_mat__) {
            for (size_t m_mat__ = 0; m_mat__ < w_m_mat_lim__; ++m_mat__) {
                w(m_mat__,n_mat__) = vals_r__[pos__++];
            }
        }
        context__.validate_dims("data initialization", "lag", "matrix_d", context__.to_vec(T,T));
        stan::math::validate_non_negative_index("lag", "T", T);
        stan::math::validate_non_negative_index("lag", "T", T);
        lag = matrix_d(T,T);
        vals_r__ = context__.vals_r("lag");
        pos__ = 0;
        size_t lag_m_mat_lim__ = T;
        size_t lag_n_mat_lim__ = T;
        for (size_t n_mat__ = 0; n_mat__ < lag_n_mat_lim__; ++n_mat__) {
            for (size_t m_mat__ = 0; m_mat__ < lag_m_mat_lim__; ++m_mat__) {
                lag(m_mat__,n_mat__) = vals_r__[pos__++];
            }
        }
        context__.validate_dims("data initialization", "N_p", "int", context__.to_vec());
        N_p = int(0);
        vals_i__ = context__.vals_i("N_p");
        pos__ = 0;
        N_p = vals_i__[pos__++];
        context__.validate_dims("data initialization", "P", "double", context__.to_vec());
        P = double(0);
        vals_r__ = context__.vals_r("P");
        pos__ = 0;
        P = vals_r__[pos__++];

        // validate data
        try { 
            check_greater_or_equal(function__,K,0,"K");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of K: ") + std::string(e.what()));
        };
        try { 
            check_greater_or_equal(function__,N,0,"N");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of N: ") + std::string(e.what()));
        };
        try { 
            check_greater_or_equal(function__,T,0,"T");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of T: ") + std::string(e.what()));
        };
        try { 
            check_greater_or_equal(function__,N_knots,0,"N_knots");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of N_knots: ") + std::string(e.what()));
        };
        try { 
            check_greater_or_equal(function__,N_cores,0,"N_cores");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of N_cores: ") + std::string(e.what()));
        };
        try { 
            check_greater_or_equal(function__,rho,0,"rho");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of rho: ") + std::string(e.what()));
        };
        try { 
            check_greater_or_equal(function__,eta,0,"eta");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of eta: ") + std::string(e.what()));
        };
        try { 
            check_greater_or_equal(function__,gamma,0,"gamma");
            check_less_or_equal(function__,gamma,1,"gamma");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of gamma: ") + std::string(e.what()));
        };
        try { 
            check_greater_or_equal(function__,psi,0,"psi");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of psi: ") + std::string(e.what()));
        };
        try { 
            check_greater_or_equal(function__,phi,0,"phi");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of phi: ") + std::string(e.what()));
        };
        try { 
            check_greater_or_equal(function__,N_p,0,"N_p");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of N_p: ") + std::string(e.what()));
        };
        W = int(0);
        stan::math::validate_non_negative_index("eta2", "K-1", K-1);
        eta2 = vector_d(K-1);
        stan::math::validate_non_negative_index("zeros", "N_knots", N_knots);
        zeros = vector_d(N_knots);
        stan::math::validate_non_negative_index("Eye_knots", "N_knots", N_knots);
        stan::math::validate_non_negative_index("Eye_knots", "N_knots", N_knots);
        Eye_knots = matrix_d(N_knots,N_knots);
        stan::math::validate_non_negative_index("C_s", "(K - 1)", (K - 1));
        stan::math::validate_non_negative_index("C_s", "N_knots", N_knots);
        stan::math::validate_non_negative_index("C_s", "N_knots", N_knots);
        C_s = std::vector<matrix_d>((K - 1),matrix_d(N_knots,N_knots));
        stan::math::validate_non_negative_index("C_s_L", "(K - 1)", (K - 1));
        stan::math::validate_non_negative_index("C_s_L", "N_knots", N_knots);
        stan::math::validate_non_negative_index("C_s_L", "N_knots", N_knots);
        C_s_L = std::vector<matrix_d>((K - 1),matrix_d(N_knots,N_knots));
        stan::math::validate_non_negative_index("C_s_inv", "(K - 1)", (K - 1));
        stan::math::validate_non_negative_index("C_s_inv", "N_knots", N_knots);
        stan::math::validate_non_negative_index("C_s_inv", "N_knots", N_knots);
        C_s_inv = std::vector<matrix_d>((K - 1),matrix_d(N_knots,N_knots));
        stan::math::validate_non_negative_index("c_s", "(K - 1)", (K - 1));
        stan::math::validate_non_negative_index("c_s", "N", N);
        stan::math::validate_non_negative_index("c_s", "N_knots", N_knots);
        c_s = std::vector<matrix_d>((K - 1),matrix_d(N,N_knots));
        stan::math::validate_non_negative_index("M", "(N * T)", (N * T));
        stan::math::validate_non_negative_index("M", "(N * T)", (N * T));
        M = matrix_d((N * T),(N * T));


	//vector<double> eta2(K-1);

	//for (int k=0; k<W; ++k) {
	 // eta2[k] = eta[k] * eta[k];
	//}
        stan::math::assign(W, (K - 1));
        for (int k = 1; k <= W; ++k) {
          stan::math::assign(get_base1_lhs(eta2,k,"eta2",1), (get_base1(eta,k,"eta",1) * get_base1(eta,k,"eta",1)));
        }
        for (int i = 1; i <= N_knots; ++i) {
            stan::math::assign(get_base1_lhs(zeros,i,"zeros",1), 0.0);
        }
        for (int i = 1; i <= N_knots; ++i) {
            for (int j = 1; j <= N_knots; ++j) {
                if (as_bool(logical_eq(i,j))) {
                    stan::math::assign(get_base1_lhs(Eye_knots,i,j,"Eye_knots",1), 1.0);
                } else {
                    stan::math::assign(get_base1_lhs(Eye_knots,i,j,"Eye_knots",1), 0);
                }
            }
        }
        for (int j = 1; j <= (N * T); ++j) {
            for (int i = 1; i <= (N * T); ++i) {
                if (as_bool(logical_eq(i,j))) {
                    stan::math::assign(get_base1_lhs(M,i,j,"M",1), (1.0 - P));
                } else {
                    stan::math::assign(get_base1_lhs(M,i,j,"M",1), -(P));
                }
            }
        }

      llt_s.resize(W);
      c_s.resize(W);
      C_s_L.resize(W);
      C_s_inv.resize(W);
      
      for (int k = 0; k < W; ++k) {
	C_s[k] = exp(-d_knots.array() / rho[k]).matrix();

	llt_s[k]      = C_s[k].llt();
	C_s_L[k]      = llt_s[k].matrixL();
	C_s_inv[k]    = llt_s[k].solve(Eigen::MatrixXd::Identity(N_knots,N_knots));
	c_s[k]        = exp(- 1.0/rho[k] * d_inter.array()).matrix();

      }

      // qr   = ones.householderQr();
      // fatQ = qr.householderQ();
      // Q    = fatQ.col(0); // thin Q, as returned in R
      // QT = Q.transpose();

        // validate transformed data

        // set parameter ranges
        num_params_r__ = 0U;
        param_ranges_i__.clear();
        ++num_params_r__;
        num_params_r__ += W;
        num_params_r__ += W;
        num_params_r__ += W;
        num_params_r__ += T * W;
        num_params_r__ += N_knots * W;
        num_params_r__ += N_knots * (W * (T - 1));
        num_params_r__ += (N * T) * W;
    }

    ~pred_model() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;


        if (!(context__.contains_r("ksi")))
            throw std::runtime_error("variable ksi missing");
        vals_r__ = context__.vals_r("ksi");
        pos__ = 0U;
        context__.validate_dims("initialization", "ksi", "double", context__.to_vec());
        double ksi(0);
        ksi = vals_r__[pos__++];
        try { writer__.scalar_lub_unconstrain(ksi_lb,ksi_ub,ksi); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable ksi: ") + e.what()); }

        if (!(context__.contains_r("sigma")))
            throw std::runtime_error("variable sigma missing");
        vals_r__ = context__.vals_r("sigma");
        pos__ = 0U;
        context__.validate_dims("initialization", "sigma", "vector_d", context__.to_vec(W));
        vector_d sigma(W);
        for (int j1__ = 0U; j1__ < W; ++j1__)
            sigma(j1__) = vals_r__[pos__++];
        try { writer__.vector_lub_unconstrain(sigma_lb,sigma_ub,sigma); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable sigma: ") + e.what()); }

        if (!(context__.contains_r("lambda")))
            throw std::runtime_error("variable lambda missing");
        vals_r__ = context__.vals_r("lambda");
        pos__ = 0U;
        context__.validate_dims("initialization", "lambda", "vector_d", context__.to_vec(W));
        vector_d lambda(W);
        for (int j1__ = 0U; j1__ < W; ++j1__)
            lambda(j1__) = vals_r__[pos__++];
        try { writer__.vector_lub_unconstrain(lambda_lb,lambda_ub,lambda); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable lambda: ") + e.what()); }

        if (!(context__.contains_r("mu")))
            throw std::runtime_error("variable mu missing");
        vals_r__ = context__.vals_r("mu");
        pos__ = 0U;
        context__.validate_dims("initialization", "mu", "vector_d", context__.to_vec(W));
        vector_d mu(W);
        for (int j1__ = 0U; j1__ < W; ++j1__)
            mu(j1__) = vals_r__[pos__++];
        try { writer__.vector_unconstrain(mu); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable mu: ") + e.what()); }

        if (!(context__.contains_r("mu_t")))
            throw std::runtime_error("variable mu_t missing");
        vals_r__ = context__.vals_r("mu_t");
        pos__ = 0U;
        context__.validate_dims("initialization", "mu_t", "vector_d", context__.to_vec(W,T));
        std::vector<vector_d> mu_t(W,vector_d(T));
        for (int j1__ = 0U; j1__ < T; ++j1__)
            for (int i0__ = 0U; i0__ < W; ++i0__)
                mu_t[i0__](j1__) = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < W; ++i0__)
            try { writer__.vector_unconstrain(mu_t[i0__]); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable mu_t: ") + e.what()); }

        if (!(context__.contains_r("alpha_s")))
            throw std::runtime_error("variable alpha_s missing");
        vals_r__ = context__.vals_r("alpha_s");
        pos__ = 0U;
        context__.validate_dims("initialization", "alpha_s", "vector_d", context__.to_vec(W,N_knots));
        std::vector<vector_d> alpha_s(W,vector_d(N_knots));
        for (int j1__ = 0U; j1__ < N_knots; ++j1__)
            for (int i0__ = 0U; i0__ < W; ++i0__)
                alpha_s[i0__](j1__) = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < W; ++i0__)
            try { writer__.vector_unconstrain(alpha_s[i0__]); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable alpha_s: ") + e.what()); }

        if (!(context__.contains_r("alpha_t")))
            throw std::runtime_error("variable alpha_t missing");
        vals_r__ = context__.vals_r("alpha_t");
        pos__ = 0U;
        context__.validate_dims("initialization", "alpha_t", "vector_d", context__.to_vec((W * (T - 1)),N_knots));
        std::vector<vector_d> alpha_t((W * (T - 1)),vector_d(N_knots));
        for (int j1__ = 0U; j1__ < N_knots; ++j1__)
            for (int i0__ = 0U; i0__ < (W * (T - 1)); ++i0__)
                alpha_t[i0__](j1__) = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < (W * (T - 1)); ++i0__)
            try { writer__.vector_unconstrain(alpha_t[i0__]); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable alpha_t: ") + e.what()); }

        if (!(context__.contains_r("g")))
            throw std::runtime_error("variable g missing");
        vals_r__ = context__.vals_r("g");
        pos__ = 0U;
        context__.validate_dims("initialization", "g", "vector_d", context__.to_vec(W,(N * T)));
        std::vector<vector_d> g(W,vector_d((N * T)));
        for (int j1__ = 0U; j1__ < (N * T); ++j1__)
            for (int i0__ = 0U; i0__ < W; ++i0__)
                g[i0__](j1__) = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < W; ++i0__)
            try { writer__.vector_unconstrain(g[i0__]); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable g: ") + e.what()); }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////


    inline double lub_transform(const double x, const double lb, const double ub,
				double &lja, double &ja, double &dj) const
    {
      double inv_logit_x;
      if (x > 0) {
        double exp_minus_x = exp(-x);
	double exp_minus_x_p1 = exp_minus_x + 1.0;
        inv_logit_x = 1.0 / (1.0 + exp_minus_x);
        lja = log(ub - lb) - x - 2 * log1p(exp_minus_x);
	ja  = (ub - lb) * exp_minus_x / (exp_minus_x_p1 * exp_minus_x_p1);
	dj  = -1.0 + 2.0 * exp_minus_x / (1 + exp_minus_x);
        if ((x < std::numeric_limits<double>::infinity())
            && (inv_logit_x == 1))
	  inv_logit_x = 1 - 1e-15;
      } else {
        double exp_x = exp(x);
	double exp_x_p1 = exp_x + 1.0;
        inv_logit_x = 1.0 - 1.0 / (1.0 + exp_x);
        lja = log(ub - lb) + x - 2 * log1p(exp_x);
	ja  = (ub - lb) * exp_x / (exp_x_p1 * exp_x_p1);
	dj  = -1.0 + 2.0 / (exp_x + 1.0);
        if ((x > -std::numeric_limits<double>::infinity())
            && (inv_logit_x== 0))
	  inv_logit_x = 1e-15;
      }
      return lb + (ub - lb) * inv_logit_x;
    }

    double normal_log_double(const vector_d& y, const double mu, const double sigma) const {
      double lp = 0.0;
      double inv_sigma = 1.0/sigma;
      double log_sigma = log(sigma);

      int size_y = y.size();

      for (int n = 0; n < size_y; n++) {
	const double y_minus_mu_over_sigma = (y[n] - mu) * inv_sigma;
	const double y_minus_mu_over_sigma_squared = y_minus_mu_over_sigma * y_minus_mu_over_sigma;
	lp -= 0.5 * y_minus_mu_over_sigma_squared;
      }

      return lp;
    }

    double normal_log_double(const double y, const double mu, const double sigma) const { 
      double lp = 0.0;
      double inv_sigma = 1.0/sigma;
      double log_sigma = log(sigma);

      const double y_minus_mu_over_sigma = (y - mu) * inv_sigma;
      const double y_minus_mu_over_sigma_squared = y_minus_mu_over_sigma * y_minus_mu_over_sigma;

      lp -= 0.5 * y_minus_mu_over_sigma_squared;
      lp -= log_sigma;

      return lp;
    }

    double multi_normal_cholesky_log_double(const vector_d& y,
					    const vector_d& mu,
					    const matrix_d& L) const {

      double lp = 0.0;

      vector_d L_log_diag = L.diagonal().array().log().matrix();
      lp -= sum(L_log_diag);

      vector_d y_minus_mu = y.array() - mu.array();
      vector_d half = mdivide_left_tri_low(L, y_minus_mu);
      lp -= 0.5 * dot_self(half);

      return lp;
    }


    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////


    template <bool propto, bool jacobian>
    double log_prob_grad(Eigen::VectorXd& params_r,
                         Eigen::VectorXd& gradient,
                         std::ostream* pstream = 0) const {

      double lp = 0.0;

      //
      // unpack model parameters and constrain
      //

      vector<double> vec_params_r;
      vector<int>    vec_params_i;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      stan::io::reader<double> in(vec_params_r, vec_params_i);

      double ksi, ksi_lja, ksi_ja, ksi_dj;
      ksi = lub_transform(in.scalar(), ksi_lb, ksi_ub, ksi_lja, ksi_ja, ksi_dj);

      if (jacobian)
      	lp += ksi_lja;

      vector_d sigma(W), sigma_lja(W), sigma_ja(W), sigma_dj(W);
      for (int i=0; i<W; ++i)
      	sigma[i] = lub_transform(in.scalar(), sigma_lb, sigma_ub, sigma_lja[i], sigma_ja[i], sigma_dj[i]);

      if (jacobian)
      	for (int i=0; i<W; ++i)
      	  lp += sigma_lja[i];

      vector_d lambda(W), lambda_lja(W), lambda_ja(W), lambda_dj(W);
      	for (int i=0; i<W; ++i)
      	  lambda[i] = lub_transform(in.scalar(), lambda_lb, lambda_ub, lambda_lja[i], lambda_ja[i], lambda_dj[i]);

      if (jacobian)
      	for (int i=0; i<W; ++i)
      	  lp += lambda_lja[i];

      vector_d  mu = in.vector_constrain(W,lp);

      vector<vector_d> mu_t(W);
      for (int k = 0; k < W; ++k) {
      	if (jacobian)
      	  mu_t[k] = in.vector_constrain(T,lp);
      	else
      	  mu_t[k] = in.vector_constrain(T);
      }

      vector<vector_d> alpha_s(W);
      for (int k = 0; k<W; ++k) {
      	if (jacobian)
      	  alpha_s[k] = in.vector_constrain(N_knots,lp);
      	else
      	  alpha_s[k] = in.vector_constrain(N_knots);
      }

      vector<vector_d> alpha_t(W*(T-1));
      for (int k = 0; k<(W*(T-1)); ++k) {
      	if (jacobian)
      	  alpha_t[k] = in.vector_constrain(N_knots,lp);
      	else
      	  alpha_t[k] = in.vector_constrain(N_knots);
      }

      vector<vector_d> g(W);
      for (int k = 0; k < W; ++k) {
      	if (jacobian)
      	  g[k] = in.vector_constrain((N * T),lp);
      	else
      	  g[k] = in.vector_constrain((N * T));
      }

      //
      // compute log probability
      //

      // priors
      lp += normal_log_double(mu, 0, 2);

      for (int k; k<W; ++k){
      	lp += normal_log_double(mu_t[k][1], 0.0, 20.0); 
      }


      std::cout << "LP after mu prior : " << lp << std::endl;

      vector<double> lambda_inv(W);
      vector<double> sigma2(W);

      vector<matrix_d> Q_s(W);
      vector<matrix_d> q_s(W);
      LLT<matrix_d> llt_Qs;
      vector<matrix_d> Q_s_L(W);
      vector<matrix_d> Q_s_inv(W);

      for (int k=0; k<W; ++k) {
      	sigma2[k] = sigma[k] * sigma[k];
      	lambda_inv[k] = 1.0 / lambda[k];
      }

      for (int k=0; k<W; ++k) {
      	Q_s[k] = (- lambda_inv[k] * d_knots.array()).exp().matrix();
      	q_s[k] = (- lambda_inv[k] * d_inter.array()).exp().matrix();
      	llt_Qs = Q_s[k].llt(); 
      	Q_s_L[k] = llt_Qs.matrixL();
      	Q_s_inv[k] = llt_Qs.solve(Eigen::MatrixXd::Identity(N_knots,N_knots));
      }


      double cvar;
      double qvar;
      
      vector_d Halpha_s;
      vector_d Halpha_t;
      vector<vector_d> qQinv_alpha(T-1);

      row_vector_d c_i;
      row_vector_d q_i;
      vector_d sqrtvar;

      vector<vector_d> mu_g(W);

      matrix_d         exp_g(N*T, W);

      for (int k; k<W; ++k){
      	lp += multi_normal_cholesky_log_double(alpha_s[k], zeros, eta2[k] * C_s_L[k]);

      	Halpha_s = c_s[k] * ( C_s_inv[k] * alpha_s[k] ); 
	
      	for (int i; i<N; ++i){
      	  for (int t; t<T; ++t){
      	    mu_g[k][i * T + t] = Halpha_s[i];
      	  }
      	}

      	lp += multi_normal_cholesky_log_double(alpha_t[k * (T-1) + 1], zeros, sigma2[k] * Q_s_L[k]);
      	for (int t=1; t<(T-1); ++t){
      	  lp += multi_normal_cholesky_log_double(alpha_t[k * (T-1) + t], zeros, sigma2[k] * Q_s_L[k]);
      	  qQinv_alpha[t] = q_s[k] * Q_s_inv[k] * alpha_t[k * (T-1) + t];
      	}
      	// time-varying mean
      	for (int t=1; t<T; ++t){
      	  lp += normal_log_double(mu_t[k][t], mu_t[k][t-1], ksi);
      	}    

      	for (int i; i<N; ++i){
      	  c_i  = c_s[k].row(i);
      	  cvar = eta2[k] * c_i * C_s_inv[k] * c_i.transpose();

      	  q_i  = q_s[k].row(i);
      	  qvar = sigma2[k] * q_i * Q_s_inv[k] * q_i.transpose();
      
      	  mu_g[k][i * T + 1] = mu[k] + mu_t[k][1] + mu_g[k][i * T + 1];
      	  sqrtvar[i * T + 1] = sqrt(eta2[k] - cvar);

      	  for (int t=1; t<T; ++i){
      	    mu_g[k][i * T + t] = mu[k] + mu_t[k][t] + mu_g[k][i * T + t] +  qQinv_alpha[t-1][i];
      	    sqrtvar[i * T + t] = sqrt(eta2[k] - cvar + sigma2[k] - qvar);	
      	  }
      	}
    
      	for (int i; i<N*T; ++i){
      	  lp += normal_log_double(g[k][i], mu_g[k][i], sqrtvar[i]);
      	}
	
    
      // //  // std::cout << "LP after alpha[k]: " <<  multi_normal_cholesky_log_double(alpha[k], zeros, C_star_L[k]) << std::endl;

      	exp_g.col(k) = exp(g[k]);
      }

      // // // for (int k=0; k<W; ++k)
      // // // 	lp += lp_thrd[k];

      // // std::cout << "LP after alpha and g: " << lp << std::endl;

      // // // #ifndef NDEBUG
      // // // if (!exp_g.allFinite()) throw std::runtime_error("exp_g is not finite");
      // // // #endif

      vector_d sum_exp_g = exp_g.rowwise().sum();
      matrix_d r(N*T,K);

      for (int k = 0; k < W; ++k)
      	for (int i = 0; i < N*T; ++i)
      	  r(i,k) = exp_g(i,k) / (1. + sum_exp_g(i));

      for (int i = 0; i < N*T; ++i)
      	r(i,W) = 1. / (1. + sum_exp_g(i));

      matrix_d r_new(N_cores * T, K);
      matrix_d out_sum(N_cores*T,K);
      vector_d sum_w(N*T);
      int      idx_core;

      out_sum.fill(0);

      // XXX: omp??
      for (int i = 0; i < N_cores; ++i) {
      	for (int t = 0; t < T; ++t) {
      	  idx_core = (idx_cores[i] - 1) * T + t;
      	  r_new.row(i * T + t) = gamma * r.row(idx_core);

      	  for (int j = 0; j < N; ++j) {
      	    if (d(idx_cores[i]-1,j) > 0) {
      	      out_sum.row(i * T + t) += w(i,j) * r.row(j * T + t);
      	    }
      	  }

      	  sum_w[i*T+t] = out_sum.row(i * T + t).sum();
      	  r_new.row(i * T + t) += out_sum.row(i *T + t) * (1 - gamma) / sum_w[i*T+t];
      	}
      }
      std::cout<<"Made it here! " << std::endl;      
      // adding these lines causes segfault...
      vector<int> N_grains(N_cores*T);
      vector_d    A(N_cores*T);
      matrix_d    kappa(N_cores*T,K);
      vector<int> y_row_sum(N_cores*T);

      // // XXX: omp?
      // for (int i = 0; i < N_cores * T; ++i) {
      // 	y_row_sum[i] = 0.0;
      // 	for (int k = 0; k < K; ++k)
      // 	  y_row_sum[i] += y[i][k];

      // 	if (y_row_sum[i] > 0) {
      // 	  for (int k = 0; k < K; ++k){
      // 	    kappa(i,k) = phi(k) * r_new(i,k);
      // 	  }

      // 	  A[i] = kappa.row(i).sum();
      // 	  N_grains[i] = y_row_sum[i];

      // 	  lp += lgamma(N_grains[i] + 1.0) + lgamma(A[i]) - lgamma(N_grains[i] + A[i]);
      // 	  for (int k = 0; k < K; ++k) {
      // 	    lp += -lgamma(y[i][k] + 1) + lgamma(y[i][k] + kappa(i,k)) - lgamma(kappa(i,k));
      // 	  }
      // 	}

      // }

      gradient.fill(0.0);

      std::cout<<"LP -----> " << lp << std::endl;

      return lp;
    }

    template <bool propto, bool jacobian>
    double log_prob(vector<double>& params_r,
    		    vector<int>& params_i,
    		    std::ostream* pstream = 0) const {

      Eigen::VectorXd evec_params_r(params_r.size());
      Eigen::VectorXd evec_gradient(params_r.size());

      for (int i=0; i<params_r.size(); i++) evec_params_r[i] = params_r[i];
      double lp = log_prob_grad<propto, jacobian>(evec_params_r, evec_gradient, pstream);
      
      return lp;
    }

    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

      throw "log_prob called";

    }

    template <bool propto, bool jacobian, typename T>
    T log_prob(Eigen::Matrix<T,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {

      throw "log_prob called";

    }


    // template <bool propto__, bool jacobian__, typename T__>
    // T__ log_prob(vector<T__>& params_r__,
    //              vector<int>& params_i__,
    //              std::ostream* pstream__ = 0) const {

    //     T__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    //     (void) DUMMY_VAR__;  // suppress unused var warning

    //     T__ lp__(0.0);
    //     stan::math::accumulator<T__> lp_accum__;

    //     // model parameters
    //     stan::io::reader<T__> in__(params_r__,params_i__);

    //     T__ ksi;
    //     (void) ksi;   // dummy to suppress unused var warning
    //     if (jacobian__)
    //         ksi = in__.scalar_lub_constrain(ksi_lb,ksi_ub,lp__);
    //     else
    //         ksi = in__.scalar_lub_constrain(ksi_lb,ksi_ub);

    //     Eigen::Matrix<T__,Eigen::Dynamic,1>  sigma;
    //     (void) sigma;   // dummy to suppress unused var warning
    //     if (jacobian__)
    //         sigma = in__.vector_lub_constrain(sigma_lb,sigma_ub,W,lp__);
    //     else
    //         sigma = in__.vector_lub_constrain(sigma_lb,sigma_ub,W);

    //     Eigen::Matrix<T__,Eigen::Dynamic,1>  lambda;
    //     (void) lambda;   // dummy to suppress unused var warning
    //     if (jacobian__)
    //         lambda = in__.vector_lb_constrain(0,W,lp__);
    //     else
    //         lambda = in__.vector_lb_constrain(0,W);

    //     Eigen::Matrix<T__,Eigen::Dynamic,1>  mu;
    //     (void) mu;   // dummy to suppress unused var warning
    //     if (jacobian__)
    //         mu = in__.vector_constrain(W,lp__);
    //     else
    //         mu = in__.vector_constrain(W);

    //     vector<Eigen::Matrix<T__,Eigen::Dynamic,1> > mu_t;
    //     size_t dim_mu_t_0__ = W;
    //     mu_t.reserve(dim_mu_t_0__);
    //     for (size_t k_0__ = 0; k_0__ < dim_mu_t_0__; ++k_0__) {
    //         if (jacobian__)
    //             mu_t.push_back(in__.vector_constrain(T,lp__));
    //         else
    //             mu_t.push_back(in__.vector_constrain(T));
    //     }

    //     vector<Eigen::Matrix<T__,Eigen::Dynamic,1> > alpha_s;
    //     size_t dim_alpha_s_0__ = W;
    //     alpha_s.reserve(dim_alpha_s_0__);
    //     for (size_t k_0__ = 0; k_0__ < dim_alpha_s_0__; ++k_0__) {
    //         if (jacobian__)
    //             alpha_s.push_back(in__.vector_constrain(N_knots,lp__));
    //         else
    //             alpha_s.push_back(in__.vector_constrain(N_knots));
    //     }

    //     vector<Eigen::Matrix<T__,Eigen::Dynamic,1> > alpha_t;
    //     size_t dim_alpha_t_0__ = (W * (T - 1));
    //     alpha_t.reserve(dim_alpha_t_0__);
    //     for (size_t k_0__ = 0; k_0__ < dim_alpha_t_0__; ++k_0__) {
    //         if (jacobian__)
    //             alpha_t.push_back(in__.vector_constrain(N_knots,lp__));
    //         else
    //             alpha_t.push_back(in__.vector_constrain(N_knots));
    //     }

    //     vector<Eigen::Matrix<T__,Eigen::Dynamic,1> > g;
    //     size_t dim_g_0__ = W;
    //     g.reserve(dim_g_0__);
    //     for (size_t k_0__ = 0; k_0__ < dim_g_0__; ++k_0__) {
    //         if (jacobian__)
    //             g.push_back(in__.vector_constrain((N * T),lp__));
    //         else
    //             g.push_back(in__.vector_constrain((N * T)));
    //     }


    //     // transformed parameters

    //     // initialized transformed params to avoid seg fault on val access
        

    //     // validate transformed parameters

    //     const char* function__ = "validate transformed params %1%";
    //     (void) function__; // dummy to suppress unused var warning
    //     // model body
    //     {
    //         vector<Eigen::Matrix<T__,Eigen::Dynamic,1> > mu_g(W, (Eigen::Matrix<T__,Eigen::Dynamic,1> ((N * T))));
    //         stan::math::fill(mu_g,DUMMY_VAR__);
    //         Eigen::Matrix<T__,Eigen::Dynamic,1>  sum_exp_g((N * T));
    //         (void) sum_exp_g;   // dummy to suppress unused var warning
    //         stan::math::fill(sum_exp_g,DUMMY_VAR__);
    //         vector<Eigen::Matrix<T__,Eigen::Dynamic,1> > r((N * T), (Eigen::Matrix<T__,Eigen::Dynamic,1> (K)));
    //         stan::math::fill(r,DUMMY_VAR__);
    //         vector<Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> > q_s(W, (Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> (N,N_knots)));
    //         stan::math::fill(q_s,DUMMY_VAR__);
    //         vector<Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> > Q_s(W, (Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> (N_knots,N_knots)));
    //         stan::math::fill(Q_s,DUMMY_VAR__);
    //         vector<Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> > Q_s_L(W, (Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> (N_knots,N_knots)));
    //         stan::math::fill(Q_s_L,DUMMY_VAR__);
    //         vector<Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> > Q_s_inv(W, (Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> (N_knots,N_knots)));
    //         stan::math::fill(Q_s_inv,DUMMY_VAR__);
    //         Eigen::Matrix<T__,Eigen::Dynamic,1>  sigma2(W);
    //         (void) sigma2;   // dummy to suppress unused var warning
    //         stan::math::fill(sigma2,DUMMY_VAR__);
    //         T__ cvar;
    //         (void) cvar;   // dummy to suppress unused var warning
    //         T__ qvar;
    //         (void) qvar;   // dummy to suppress unused var warning
    //         Eigen::Matrix<T__,Eigen::Dynamic,1>  Halpha_s(N);
    //         (void) Halpha_s;   // dummy to suppress unused var warning
    //         stan::math::fill(Halpha_s,DUMMY_VAR__);
    //         Eigen::Matrix<T__,Eigen::Dynamic,1>  Halpha_t(N);
    //         (void) Halpha_t;   // dummy to suppress unused var warning
    //         stan::math::fill(Halpha_t,DUMMY_VAR__);
    //         vector<Eigen::Matrix<T__,Eigen::Dynamic,1> > qQinv_alpha((T - 1), (Eigen::Matrix<T__,Eigen::Dynamic,1> (N)));
    //         stan::math::fill(qQinv_alpha,DUMMY_VAR__);
    //         Eigen::Matrix<T__,1,Eigen::Dynamic>  c_i(N_knots);
    //         (void) c_i;   // dummy to suppress unused var warning
    //         stan::math::fill(c_i,DUMMY_VAR__);
    //         Eigen::Matrix<T__,1,Eigen::Dynamic>  q_i(N_knots);
    //         (void) q_i;   // dummy to suppress unused var warning
    //         stan::math::fill(q_i,DUMMY_VAR__);
    //         Eigen::Matrix<T__,Eigen::Dynamic,1>  sqrtvar((N * T));
    //         (void) sqrtvar;   // dummy to suppress unused var warning
    //         stan::math::fill(sqrtvar,DUMMY_VAR__);
    //         stan::math::initialize(mu_g, DUMMY_VAR__);
    //         stan::math::initialize(sum_exp_g, DUMMY_VAR__);
    //         stan::math::initialize(r, DUMMY_VAR__);
    //         stan::math::initialize(q_s, DUMMY_VAR__);
    //         stan::math::initialize(Q_s, DUMMY_VAR__);
    //         stan::math::initialize(Q_s_L, DUMMY_VAR__);
    //         stan::math::initialize(Q_s_inv, DUMMY_VAR__);
    //         stan::math::initialize(sigma2, DUMMY_VAR__);
    //         stan::math::initialize(cvar, DUMMY_VAR__);
    //         stan::math::initialize(qvar, DUMMY_VAR__);
    //         stan::math::initialize(Halpha_s, DUMMY_VAR__);
    //         stan::math::initialize(Halpha_t, DUMMY_VAR__);
    //         stan::math::initialize(qQinv_alpha, DUMMY_VAR__);
    //         stan::math::initialize(c_i, DUMMY_VAR__);
    //         stan::math::initialize(q_i, DUMMY_VAR__);
    //         stan::math::initialize(sqrtvar, DUMMY_VAR__);
    //         lp_accum__.add(normal_log<propto__>(mu, 0, 20));
    //         for (int k = 1; k <= W; ++k) {
    //             lp_accum__.add(normal_log<propto__>(get_base1(get_base1(mu_t,k,"mu_t",1),1,"mu_t",2), 0.0, 20.0));
    //         }
    //         for (int k = 1; k <= W; ++k) {
    //             stan::math::assign(get_base1_lhs(sigma2,k,"sigma2",1), (get_base1(sigma,k,"sigma",1) * get_base1(sigma,k,"sigma",1)));
    //         }
    //         for (int k = 1; k <= W; ++k) {
    //             stan::math::assign(get_base1_lhs(q_s,k,"q_s",1), exp(divide(minus(d_inter),get_base1(lambda,k,"lambda",1))));
    //             stan::math::assign(get_base1_lhs(Q_s,k,"Q_s",1), exp(divide(minus(d_knots),get_base1(lambda,k,"lambda",1))));
    //             stan::math::assign(get_base1_lhs(Q_s_L,k,"Q_s_L",1), cholesky_decompose(get_base1(Q_s,k,"Q_s",1)));
    //             stan::math::assign(get_base1_lhs(Q_s_inv,k,"Q_s_inv",1), transpose(mdivide_right_tri_low(transpose(mdivide_left_tri_low(get_base1(Q_s_L,k,"Q_s_L",1),Eye_knots)),get_base1(Q_s_L,k,"Q_s_L",1))));
    //         }
    //         for (int k = 1; k <= W; ++k) {
    //             lp_accum__.add(multi_normal_prec_log<propto__>(get_base1(alpha_s,k,"alpha_s",1), zeros, multiply((1 / get_base1(eta2,k,"eta2",1)),get_base1(C_s_inv,k,"C_s_inv",1))));
    //             stan::math::assign(Halpha_s, multiply(get_base1(c_s,k,"c_s",1),multiply(get_base1(C_s_inv,k,"C_s_inv",1),get_base1(alpha_s,k,"alpha_s",1))));
    //             for (int i = 1; i <= N; ++i) {
    //                 for (int t = 1; t <= T; ++t) {
    //                     stan::math::assign(get_base1_lhs(get_base1_lhs(mu_g,k,"mu_g",1),(((i - 1) * T) + t),"mu_g",2), get_base1(Halpha_s,i,"Halpha_s",1));
    //                 }
    //             }
    //             stan::math::assign(get_base1_lhs(mu_g,k,"mu_g",1), get_base1(mu_g,k,"mu_g",1));
    //             if (pstream__) {
    //                 stan_print(pstream__,"Here");
    //                 *pstream__ << std::endl;
    //             }
    //             lp_accum__.add(multi_normal_prec_log<propto__>(get_base1(alpha_t,(((k - 1) * (T - 1)) + 1),"alpha_t",1), zeros, multiply((1 / get_base1(sigma2,k,"sigma2",1)),get_base1(Q_s_inv,k,"Q_s_inv",1))));
    //             if (pstream__) {
    //                 stan_print(pstream__,"k = ");
    //                 stan_print(pstream__,k);
    //                 *pstream__ << std::endl;
    //             }
    //             for (int t = 2; t <= (T - 1); ++t) {
    //                 if (pstream__) {
    //                     stan_print(pstream__,"index = ");
    //                     stan_print(pstream__,((((k - 1) * (T - 1)) + t) - 1));
    //                     *pstream__ << std::endl;
    //                 }
    //                 lp_accum__.add(multi_normal_prec_log<propto__>(get_base1(alpha_t,(((k - 1) * (T - 1)) + t),"alpha_t",1), get_base1(alpha_t,((((k - 1) * (T - 1)) + t) - 1),"alpha_t",1), multiply((1 / get_base1(sigma2,k,"sigma2",1)),get_base1(Q_s_inv,k,"Q_s_inv",1))));
    //             }
    //             for (int t = 1; t <= (T - 1); ++t) {
    //                 stan::math::assign(get_base1_lhs(qQinv_alpha,t,"qQinv_alpha",1), multiply(multiply(get_base1(q_s,k,"q_s",1),get_base1(Q_s_inv,k,"Q_s_inv",1)),get_base1(alpha_t,(((k - 1) * (T - 1)) + t),"alpha_t",1)));
    //             }
    //             for (int i = 2; i <= T; ++i) {
    //                 lp_accum__.add(normal_log<propto__>(get_base1(get_base1(mu_t,k,"mu_t",1),i,"mu_t",2), get_base1(get_base1(mu_t,k,"mu_t",1),(i - 1),"mu_t",2), ksi));
    //             }
    //             for (int i = 1; i <= N; ++i) {
    //                 stan::math::assign(c_i, row(get_base1(c_s,k,"c_s",1),i));
    //                 stan::math::assign(cvar, multiply(multiply(multiply(get_base1(eta2,k,"eta2",1),c_i),get_base1(C_s_inv,k,"C_s_inv",1)),transpose(c_i)));
    //                 stan::math::assign(q_i, row(get_base1(q_s,k,"q_s",1),i));
    //                 stan::math::assign(qvar, multiply(multiply(multiply(get_base1(sigma2,k,"sigma2",1),q_i),get_base1(Q_s_inv,k,"Q_s_inv",1)),transpose(q_i)));
    //                 stan::math::assign(get_base1_lhs(get_base1_lhs(mu_g,k,"mu_g",1),(((i - 1) * T) + 1),"mu_g",2), ((get_base1(mu,k,"mu",1) + get_base1(get_base1(mu_t,k,"mu_t",1),1,"mu_t",2)) + get_base1(get_base1(mu_g,k,"mu_g",1),(((i - 1) * T) + 1),"mu_g",2)));
    //                 stan::math::assign(get_base1_lhs(sqrtvar,(((i - 1) * T) + 1),"sqrtvar",1), sqrt((get_base1(eta2,k,"eta2",1) - cvar)));
    //                 for (int t = 2; t <= T; ++t) {
    //                     stan::math::assign(get_base1_lhs(get_base1_lhs(mu_g,k,"mu_g",1),(((i - 1) * T) + t),"mu_g",2), (((get_base1(mu,k,"mu",1) + get_base1(get_base1(mu_t,k,"mu_t",1),t,"mu_t",2)) + get_base1(get_base1(mu_g,k,"mu_g",1),(((i - 1) * T) + t),"mu_g",2)) + get_base1(get_base1(qQinv_alpha,(t - 1),"qQinv_alpha",1),i,"qQinv_alpha",2)));
    //                     stan::math::assign(get_base1_lhs(sqrtvar,(((i - 1) * T) + t),"sqrtvar",1), sqrt((((get_base1(eta2,k,"eta2",1) - cvar) + get_base1(sigma2,k,"sigma2",1)) - qvar)));
    //                 }
    //             }
    //             for (int i = 1; i <= (N * T); ++i) {
    //                 lp_accum__.add(normal_log<propto__>(get_base1(get_base1(g,k,"g",1),i,"g",2), get_base1(get_base1(mu_g,k,"mu_g",1),i,"mu_g",2), get_base1(sqrtvar,i,"sqrtvar",1)));
    //             }
    //         }
    //         for (int i = 1; i <= (N * T); ++i) {
    //             stan::math::assign(get_base1_lhs(sum_exp_g,i,"sum_exp_g",1), 0.0);
    //             for (int k = 1; k <= W; ++k) {
    //                 stan::math::assign(get_base1_lhs(sum_exp_g,i,"sum_exp_g",1), (get_base1(sum_exp_g,i,"sum_exp_g",1) + exp(get_base1(get_base1(g,k,"g",1),i,"g",2))));
    //             }
    //         }
    //         for (int k = 1; k <= W; ++k) {
    //             for (int i = 1; i <= (N * T); ++i) {
    //                 stan::math::assign(get_base1_lhs(get_base1_lhs(r,i,"r",1),k,"r",2), (exp(get_base1(get_base1(g,k,"g",1),i,"g",2)) / (1 + get_base1(sum_exp_g,i,"sum_exp_g",1))));
    //             }
    //         }
    //         for (int i = 1; i <= (N * T); ++i) {
    //             stan::math::assign(get_base1_lhs(get_base1_lhs(r,i,"r",1),K,"r",2), (1 / (1 + get_base1(sum_exp_g,i,"sum_exp_g",1))));
    //         }
    //         {
    //             vector<Eigen::Matrix<T__,Eigen::Dynamic,1> > r_new((N_cores * T), (Eigen::Matrix<T__,Eigen::Dynamic,1> (K)));
    //             stan::math::fill(r_new,DUMMY_VAR__);
    //             Eigen::Matrix<T__,Eigen::Dynamic,1>  out_sum(K);
    //             (void) out_sum;   // dummy to suppress unused var warning
    //             stan::math::fill(out_sum,DUMMY_VAR__);
    //             T__ sum_w;
    //             (void) sum_w;   // dummy to suppress unused var warning
    //             int idx_core(0);
    //             (void) idx_core;   // dummy to suppress unused var warning
    //             stan::math::initialize(r_new, DUMMY_VAR__);
    //             stan::math::initialize(out_sum, DUMMY_VAR__);
    //             stan::math::initialize(sum_w, DUMMY_VAR__);
    //             for (int i = 1; i <= N_cores; ++i) {
    //                 for (int t = 1; t <= T; ++t) {
    //                     stan::math::assign(idx_core, (((get_base1(idx_cores,i,"idx_cores",1) - 1) * T) + t));
    //                     stan::math::assign(get_base1_lhs(r_new,(((i - 1) * T) + t),"r_new",1), multiply(gamma,get_base1(r,idx_core,"r",1)));
    //                     for (int k = 1; k <= K; ++k) {
    //                         stan::math::assign(get_base1_lhs(out_sum,k,"out_sum",1), 0);
    //                     }
    //                     stan::math::assign(sum_w, 0);
    //                     for (int j = 1; j <= N; ++j) {
    //                         if (as_bool(logical_gt(get_base1(d,get_base1(idx_cores,i,"idx_cores",1),j,"d",1),0))) {
    //                             stan::math::assign(out_sum, add(out_sum,multiply(get_base1(w,i,j,"w",1),get_base1(r,(((j - 1) * T) + t),"r",1))));
    //                         }
    //                     }
    //                     stan::math::assign(sum_w, sum(out_sum));
    //                     stan::math::assign(get_base1_lhs(r_new,(((i - 1) * T) + t),"r_new",1), add(get_base1(r_new,(((i - 1) * T) + t),"r_new",1),divide(multiply(out_sum,(1 - gamma)),sum_w)));
    //                 }
    //             }
    //             {
    //                 T__ N_grains;
    //                 (void) N_grains;   // dummy to suppress unused var warning
    //                 T__ A;
    //                 (void) A;   // dummy to suppress unused var warning
    //                 Eigen::Matrix<T__,Eigen::Dynamic,1>  kappa(K);
    //                 (void) kappa;   // dummy to suppress unused var warning
    //                 stan::math::fill(kappa,DUMMY_VAR__);
    //                 stan::math::initialize(N_grains, DUMMY_VAR__);
    //                 stan::math::initialize(A, DUMMY_VAR__);
    //                 stan::math::initialize(kappa, DUMMY_VAR__);
    //                 for (int i = 1; i <= (N_cores * T); ++i) {
    //                     if (as_bool(logical_gt(sum(get_base1(y,i,"y",1)),0))) {
    //                         stan::math::assign(kappa, elt_multiply(phi,get_base1(r_new,i,"r_new",1)));
    //                         stan::math::assign(A, sum(kappa));
    //                         stan::math::assign(N_grains, sum(get_base1(y,i,"y",1)));
    //                         lp_accum__.add(((lgamma((N_grains + 1)) + lgamma(A)) - lgamma((N_grains + A))));
    //                         for (int k = 1; k <= K; ++k) {
    //                             lp_accum__.add(((-(lgamma((get_base1(get_base1(y,i,"y",1),k,"y",2) + 1))) + lgamma((get_base1(get_base1(y,i,"y",1),k,"y",2) + get_base1(kappa,k,"kappa",1)))) - lgamma(get_base1(kappa,k,"kappa",1))));
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }

    //     lp_accum__.add(lp__);
    //     return lp_accum__.sum();

    // } // log_prob()

    // template <bool propto, bool jacobian, typename T>
    // T log_prob(Eigen::Matrix<T,Eigen::Dynamic,1>& params_r,
    //            std::ostream* pstream = 0) const {
    //   std::vector<T> vec_params_r;
    //   vec_params_r.reserve(params_r.size());
    //   for (int i = 0; i < params_r.size(); ++i)
    //     vec_params_r.push_back(params_r(i));
    //   std::vector<int> vec_params_i;
    //   return log_prob<propto,jacobian,T>(vec_params_r, vec_params_i, pstream);
    // }


    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("ksi");
        names__.push_back("sigma");
        names__.push_back("lambda");
        names__.push_back("mu");
        names__.push_back("mu_t");
        names__.push_back("alpha_s");
        names__.push_back("alpha_t");
        names__.push_back("g");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(W);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(W);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(W);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(W);
        dims__.push_back(T);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(W);
        dims__.push_back(N_knots);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back((W * (T - 1)));
        dims__.push_back(N_knots);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(W);
        dims__.push_back((N * T));
        dimss__.push_back(dims__);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        vars__.resize(0);
        stan::io::reader<double> in__(params_r__,params_i__);
        static const char* function__ = "pred_model_namespace::write_array(%1%)";
        (void) function__; // dummy call to supress warning
        // read-transform, write parameters
        double ksi = in__.scalar_lub_constrain(ksi_lb,ksi_ub);
        vector_d sigma = in__.vector_lub_constrain(sigma_lb,sigma_ub,W);
        vector_d lambda = in__.vector_lub_constrain(lambda_lb,lambda_ub,W);
        vector_d mu = in__.vector_constrain(W);
        vector<vector_d> mu_t;
        size_t dim_mu_t_0__ = W;
        for (size_t k_0__ = 0; k_0__ < dim_mu_t_0__; ++k_0__) {
            mu_t.push_back(in__.vector_constrain(T));
        }
        vector<vector_d> alpha_s;
        size_t dim_alpha_s_0__ = W;
        for (size_t k_0__ = 0; k_0__ < dim_alpha_s_0__; ++k_0__) {
            alpha_s.push_back(in__.vector_constrain(N_knots));
        }
        vector<vector_d> alpha_t;
        size_t dim_alpha_t_0__ = (W * (T - 1));
        for (size_t k_0__ = 0; k_0__ < dim_alpha_t_0__; ++k_0__) {
            alpha_t.push_back(in__.vector_constrain(N_knots));
        }
        vector<vector_d> g;
        size_t dim_g_0__ = W;
        for (size_t k_0__ = 0; k_0__ < dim_g_0__; ++k_0__) {
            g.push_back(in__.vector_constrain((N * T)));
        }
        vars__.push_back(ksi);
        for (int k_0__ = 0; k_0__ < W; ++k_0__) {
            vars__.push_back(sigma[k_0__]);
        }
        for (int k_0__ = 0; k_0__ < W; ++k_0__) {
            vars__.push_back(lambda[k_0__]);
        }
        for (int k_0__ = 0; k_0__ < W; ++k_0__) {
            vars__.push_back(mu[k_0__]);
        }
        for (int k_1__ = 0; k_1__ < T; ++k_1__) {
            for (int k_0__ = 0; k_0__ < W; ++k_0__) {
                vars__.push_back(mu_t[k_0__][k_1__]);
            }
        }
        for (int k_1__ = 0; k_1__ < N_knots; ++k_1__) {
            for (int k_0__ = 0; k_0__ < W; ++k_0__) {
                vars__.push_back(alpha_s[k_0__][k_1__]);
            }
        }
        for (int k_1__ = 0; k_1__ < N_knots; ++k_1__) {
            for (int k_0__ = 0; k_0__ < (W * (T - 1)); ++k_0__) {
                vars__.push_back(alpha_t[k_0__][k_1__]);
            }
        }
        for (int k_1__ = 0; k_1__ < (N * T); ++k_1__) {
            for (int k_0__ = 0; k_0__ < W; ++k_0__) {
                vars__.push_back(g[k_0__][k_1__]);
            }
        }

        if (!include_tparams__) return;
        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__; // dummy call to supress warning
        stan::math::accumulator<double> lp_accum__;



        // validate transformed parameters

        // write transformed parameters

        if (!include_gqs__) return;
        // declare and define generated quantities


        // validate generated quantities

        // write generated quantities
    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }


    void write_csv_header(std::ostream& o__) const {
        stan::io::csv_writer writer__(o__);
        writer__.comma();
        o__ << "ksi";
        for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
            writer__.comma();
            o__ << "sigma" << '.' << k_0__;
        }
        for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
            writer__.comma();
            o__ << "lambda" << '.' << k_0__;
        }
        for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
            writer__.comma();
            o__ << "mu" << '.' << k_0__;
        }
        for (int k_1__ = 1; k_1__ <= T; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
                writer__.comma();
                o__ << "mu_t" << '.' << k_0__ << '.' << k_1__;
            }
        }
        for (int k_1__ = 1; k_1__ <= N_knots; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
                writer__.comma();
                o__ << "alpha_s" << '.' << k_0__ << '.' << k_1__;
            }
        }
        for (int k_1__ = 1; k_1__ <= N_knots; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= (W * (T - 1)); ++k_0__) {
                writer__.comma();
                o__ << "alpha_t" << '.' << k_0__ << '.' << k_1__;
            }
        }
        for (int k_1__ = 1; k_1__ <= (N * T); ++k_1__) {
            for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
                writer__.comma();
                o__ << "g" << '.' << k_0__ << '.' << k_1__;
            }
        }
        writer__.newline();
    }

    template <typename RNG>
    void write_csv(RNG& base_rng__,
                   std::vector<double>& params_r__,
                   std::vector<int>& params_i__,
                   std::ostream& o__,
                   std::ostream* pstream__ = 0) const {
        stan::io::reader<double> in__(params_r__,params_i__);
        stan::io::csv_writer writer__(o__);
        static const char* function__ = "pred_model_namespace::write_csv(%1%)";
        (void) function__; // dummy call to supress warning
        // read-transform, write parameters
        double ksi = in__.scalar_lub_constrain(ksi_lb,ksi_ub);
        writer__.write(ksi);
        vector_d sigma = in__.vector_lub_constrain(sigma_lb,sigma_ub,W);
        writer__.write(sigma);
        vector_d lambda = in__.vector_lub_constrain(lambda_lb,lambda_ub, W);
        writer__.write(lambda);
        vector_d mu = in__.vector_constrain(W);
        writer__.write(mu);
        vector<vector_d> mu_t;
        size_t dim_mu_t_0__ = W;
        for (size_t k_0__ = 0; k_0__ < dim_mu_t_0__; ++k_0__) {
            mu_t.push_back(in__.vector_constrain(T));
            writer__.write(mu_t[k_0__]);
        }
        vector<vector_d> alpha_s;
        size_t dim_alpha_s_0__ = W;
        for (size_t k_0__ = 0; k_0__ < dim_alpha_s_0__; ++k_0__) {
            alpha_s.push_back(in__.vector_constrain(N_knots));
            writer__.write(alpha_s[k_0__]);
        }
        vector<vector_d> alpha_t;
        size_t dim_alpha_t_0__ = (W * (T - 1));
        for (size_t k_0__ = 0; k_0__ < dim_alpha_t_0__; ++k_0__) {
            alpha_t.push_back(in__.vector_constrain(N_knots));
            writer__.write(alpha_t[k_0__]);
        }
        vector<vector_d> g;
        size_t dim_g_0__ = W;
        for (size_t k_0__ = 0; k_0__ < dim_g_0__; ++k_0__) {
            g.push_back(in__.vector_constrain((N * T)));
            writer__.write(g[k_0__]);
        }

        // declare, define and validate transformed parameters
        double lp__ = 0.0;
        (void) lp__; // dummy call to supress warning
        stan::math::accumulator<double> lp_accum__;




        // write transformed parameters

        // declare and define generated quantities


        // validate generated quantities

        // write generated quantities
        writer__.newline();
    }

    template <typename RNG>
    void write_csv(RNG& base_rng,
                   Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                   std::ostream& o,
                   std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<int> params_i_vec;  // dummy
      write_csv(base_rng, params_r_vec, params_i_vec, o, pstream);
    }

    static std::string model_name() {
        return "pred_model";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "ksi";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_1__ = 1; k_1__ <= T; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu_t" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_1__ = 1; k_1__ <= N_knots; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "alpha_s" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_1__ = 1; k_1__ <= N_knots; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= (W * (T - 1)); ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "alpha_t" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_1__ = 1; k_1__ <= (N * T); ++k_1__) {
            for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "g" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "ksi";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_1__ = 1; k_1__ <= T; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu_t" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_1__ = 1; k_1__ <= N_knots; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "alpha_s" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_1__ = 1; k_1__ <= N_knots; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= (W * (T - 1)); ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "alpha_t" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_1__ = 1; k_1__ <= (N * T); ++k_1__) {
            for (int k_0__ = 1; k_0__ <= W; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "g" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
    }

}; // model

} // namespace

namespace stan {
  namespace model {

    template <bool propto, bool jacobian_adjust_transform>
    double log_prob_grad(const pred_model_namespace::pred_model& model,
                         Eigen::VectorXd& params_r,
                         Eigen::VectorXd& gradient,
                         std::ostream* msgs = 0) {
      double lp = model.template log_prob_grad<propto, jacobian_adjust_transform>(params_r, gradient, msgs);
      return lp;

    }

    template <bool propto, bool jacobian_adjust_transform>
    double log_prob_grad(const pred_model_namespace::pred_model& model,
                         std::vector<double>& params_r,
                         std::vector<int>& params_i,
                         std::vector<double>& gradient,
                         std::ostream* msgs = 0) {

      Eigen::VectorXd evec_params_r(params_r.size());
      Eigen::VectorXd evec_gradient(params_r.size());

      for (int i=0; i<params_r.size(); i++) evec_params_r[i] = params_r[i];
      double lp = model.template log_prob_grad<propto, jacobian_adjust_transform>(evec_params_r, evec_gradient, msgs);

      gradient.resize(params_r.size());
      for (int i=0; i<params_r.size(); i++) gradient[i] = evec_gradient[i];
      return lp;

    }

  }
}

#include <stan/gm/command.hpp>

int main(int argc, const char* argv[]) {
    try {
        return stan::gm::command<pred_model_namespace::pred_model>(argc,argv);
    } catch (std::exception& e) {
        std::cerr << std::endl << "Exception: " << e.what() << std::endl;
        std::cerr << "Diagnostic information: " << std::endl << boost::diagnostic_information(e) << std::endl;
        return -1;
    }
}

