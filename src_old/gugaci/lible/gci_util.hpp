#pragma once

#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <lible/configuration.hpp>
#include <lible/wave_function.hpp>

#include <ankerl/unordered_dense.h> 
#include <armadillo>

namespace lible
{
    namespace guga
    {
        using nonet = std::tuple<int, int, int, int, int, int, int, int, int>;
        using quintet = std::tuple<int, int, int, int, int>;

        using connection_list_dia = std::vector<std::tuple<int, int>>;
        using connection_list_1el = std::vector<std::tuple<int, int, int, bool>>;
        using connection_list_2el = std::vector<std::tuple<int, int, int, bool, bool>>;

        using connection_map_1el = std::map<quintet, connection_list_1el>;
        using connection_map_2el = std::map<nonet, connection_list_2el>;
        using connection_map_dia = std::map<quintet, connection_list_dia>;
        using connection_map_EpqErr = std::map<quintet,
                                               std::map<std::tuple<int, int, bool>,
                                                        std::vector<std::tuple<int, int>>>>;

        using cc_map = ankerl::unordered_dense::map<std::pair<int, int>, double>;

        using proto_1el_tuple = std::tuple<CFGProto, CFGProto>;
        using proto_2el_tuple = std::tuple<CFGProto, CFGProto, CFGProto, CFGProto>;

        using sfs_pair_t = std::pair<std::set<int>, std::set<int>>;
        using sf_pair_map_1el = std::map<quintet, sfs_pair_t>;
        using sf_pair_map_2el = std::map<nonet, sfs_pair_t>;

        using wfn_ptr = std::unique_ptr<WaveFunction>;

        struct DataFOIS
        {
            DataFOIS()
            {
                wfn = std::make_unique<WaveFunction>();
            }

            std::map<size_t, std::set<std::string>> cfgs_sfs;
            connection_map_1el connections_1el;
            connection_map_2el connections_2el;
            connection_map_EpqErr connections_EpqErr;
            connection_map_EpqErr connections_ErrEpq;
            wfn_ptr wfn;            
        };

        struct DataVar
        {
            DataVar()
            {
                wfn = std::make_unique<WaveFunction>();
            }

            std::map<size_t, std::set<std::string>> cfgs_sfs;
            connection_map_1el connections_1el;
            connection_map_2el connections_2el;
            connection_map_EpqErr connections_EpqErr;
            connection_map_EpqErr connections_ErrEpq;
            wfn_ptr wfn;
        };

        enum ExcType
        {
            DS,
            DV,
            SS,
            SV
        };

        namespace util
        {
            double A(const int &b, const int &x, const int &y);
            double C(const int &b, const int &x);
            double f(const char &d, const int &b);
            double fNS(const int &nue, const double &spin);
            int determineStepb(const char &d);

            bool determine1ElPhase(const int &p, const int &q, const std::string &cfg_right);

            bool determine2ElPhase(const int &p, const int &q, const int &r, const int &s,
                                   const std::string &conf_ri, const std::string &conf_right);

            size_t pq2DTo1D(const int &p, const int &q, const int &n_orb);

            size_t pqrs4DTo1D(const int &p, const int &q, const int &r, const int &s,
                              const int &n_orb);

            std::tuple<int, int> pq1DTo2D(const size_t &pq, const int &n_orbs);

            std::tuple<int, int, int, int> pqrs1DTo4D(const size_t &pqrs, const int &n_orbs);

            std::tuple<std::string, std::string> extractCFGandSF(const std::string &csf);

            std::string extractSF(const std::string &csf);

            std::vector<std::string> returnSFs(const std::map<int, std::string> &sf_map,
                                     const std::vector<int> &sf_idxs);

            std::map<std::string, int> returnSFMap(const std::map<std::string, int> &sf_map_in,
                                                   const std::set<std::string> &sfs);

            quintet returnCCInfo(const int &p, const int &q, const int &nue_left,
                                 const int &nue_right, const std::string &cfg_left,
                                 const std::string &cfg_right);

            std::vector<std::string> findConnectedSFs_DOMOSOMO(const int &p, const int &q,
                                                               const CFGProto &right);

            std::vector<std::string> findConnectedSFs_DOMOVirtual(const int &p, const int &q,
                                                                  const CFGProto &right);

            std::vector<std::string> findConnectedSFs_SOMOSOMO(const int &p, const int &q,
                                                               const CFGProto &right);

            std::vector<std::string> findConnectedSFs_SOMOVirtual(const int &p, const int &q,
                                                                  const CFGProto &right);

            arma::dmat calcCCDOMOSOMO(const int &prel, const int &qrel,
                                      const CFGProto &cfg_left,
                                      const CFGProto &cfg_right);

            arma::dmat calcCCDOMOVirtual(const int &prel, const int &qrel,
                                         const CFGProto &cfg_left,
                                         const CFGProto &cfg_right);

            arma::dmat calcCCSOMOSOMO(const int &prel, const int &qrel,
                                      const CFGProto &cfg_left,
                                      const CFGProto &cfg_right);

            arma::dmat calcCCSOMOVirtual(const int &prel, const int &qrel,
                                         const CFGProto &cfg_left,
                                         const CFGProto &cfg_right);

            void climbShavittGraphs_DOMOSOMO_R(const int &p, const int &q, bool flip,
                                               int b_left, int i, std::string sf,
                                               CSFTree::Node *node_right,
                                               std::set<std::string> &connected_sfs);

            void climbShavittGraphs_DOMOSOMO_L(const int &p, const int &q, bool flip,
                                               int b_left, int i, std::string sf,
                                               CSFTree::Node *node_right,
                                               std::set<std::string> &connected_sfs);

            void climbShavittGraphs_DOMOVirtual_R(const int &p, const int &q, bool flip,
                                                  int b_right, int i, std::string sf,
                                                  CSFTree::Node *node_left,
                                                  std::set<std::string> &connected_sfs);

            void climbShavittGraphs_DOMOVirtual_L(const int &p, const int &q, bool flip,
                                                  int b_right, int i, std::string sf,
                                                  CSFTree::Node *node_left,
                                                  std::set<std::string> &connected_sfs);

            void climbShavittGraphs_SOMOSOMO_R(const int &p, const int &q, bool flip,
                                               int b_left, int i, std::string sf,
                                               CSFTree::Node *node_right,
                                               std::set<std::string> &connected_sfs);

            void climbShavittGraphs_SOMOSOMO_L(const int &p, const int &q, bool flip,
                                               int b_left, int i, std::string sf,
                                               CSFTree::Node *node_right,
                                               std::set<std::string> &connected_sfs);

            void climbShavittGraphs_SOMOVirtual_R(const int &p, const int &q, bool flip,
                                                  int b_left, int i, std::string sf,
                                                  CSFTree::Node *node_right,
                                                  std::set<std::string> &connected_sfs);

            void climbShavittGraphs_SOMOVirtual_L(const int &p, const int &q, bool flip,
                                                  int b_left, int i, std::string sf,
                                                  CSFTree::Node *node_right,
                                                  std::set<std::string> &connected_sfs);

            void walkShavittGraphs_DOMOSOMO_R(const int &p, const int &q,
                                              bool flip, double cc, int i,
                                              CSFTree::Node *node_left,
                                              CSFTree::Node *node_right,
                                              arma::dmat &cc_matrix);

            void walkShavittGraphs_DOMOSOMO_L(const int &p, const int &q,
                                              bool flip, double cc, int i,
                                              CSFTree::Node *node_left,
                                              CSFTree::Node *node_right,
                                              arma::dmat &cc_matrix);

            void walkShavittGraphs_SOMOSOMO_R(const int &p, const int &q,
                                              bool flip, double cc, int i,
                                              CSFTree::Node *node_left,
                                              CSFTree::Node *node_right,
                                              arma::dmat &cc_matrix);

            void walkShavittGraphs_SOMOSOMO_L(const int &p, const int &q,
                                              bool flip, double cc, int i,
                                              CSFTree::Node *node_left,
                                              CSFTree::Node *node_right,
                                              arma::dmat &cc_matrix);

            void walkShavittGraphs_SOMOVirtual_R(const int &p, const int &q,
                                                 bool flip, double cc, int i,
                                                 CSFTree::Node *node_left,
                                                 CSFTree::Node *node_right,
                                                 arma::dmat &cc_matrix);

            void walkShavittGraphs_SOMOVirtual_L(const int &p, const int &q,
                                                 bool flip, double cc, int i,
                                                 CSFTree::Node *node_left,
                                                 CSFTree::Node *node_right,
                                                 arma::dmat &cc_matrix);

            template <typename T>
            void mergeSFPairs(T &sf_pairs_in, T &sf_pairs_out);

            template <typename T>
            void mergeConnections(const T &connections_in, T &connections_out);

            template <typename T>
            void mergeCCs(const T &ccs_in, T &ccs_out);
        }
    }
}
