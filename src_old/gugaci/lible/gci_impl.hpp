#pragma once

#include <memory>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include <lible/gci.hpp>
#include <lible/gci_settings.hpp>
#include <lible/gci_util.hpp>

#ifdef _USE_MPI_
#include <mpi.h>
#endif

namespace lible
{
	namespace guga
	{
		class GCI::Impl
		{
		public:
			Impl(const int &n_orbs, const int &n_els,
				 const int &n_roots, const int &multiplicity,
				 const vec2d &one_el_ints, const vec4d &two_el_ints,
				 const double &core_energy);
			~Impl();

			Impl(const Impl &) = delete;
			Impl &operator=(const Impl &) = delete;

			Impl(Impl &&) noexcept;
			Impl &operator=(Impl &&) noexcept;

			void run(std::vector<double> &ci_energies_out,
					 std::vector<std::vector<double>> &ci_vectors_out);

			// void runFromCFGs(const std::vector<std::string> &cfgs,
			// 				 std::vector<double> &ci_energies_out,
			// 				 std::vector<std::vector<double>> &ci_vectors_out);

			// void runFromCSFs(const std::vector<std::string> &csfs,
			// 				 std::vector<double> &ci_energies_out,
			// 				 std::vector<std::vector<double>> &ci_vectors_out);

			// void runFromCFGsFile(const std::string &cfgs_fname,
			// 					 std::vector<double> &ci_energies_out,
			// 					 std::vector<std::vector<double>> &ci_vectors_out);

			void runFromCSFsFile(const std::string &csfs_fname,
								 std::vector<double> &ci_energies_out,
								 std::vector<std::vector<double>> &ci_vectors_out);

			void calc1RDM(std::vector<double> &one_rdm_out);

			void calc2RDM(std::vector<double> &two_rdm_out);

			void calcSpin1RDM(std::vector<double> &one_srdm_out);

			std::vector<double> calcSigma(const std::vector<double> &trial);
			
			std::vector<double> calcSigma(const vec2d &aux_1el_ints,
										  const vec4d &aux_2el_ints,
										  const std::vector<double> &trial);

		private:
			/* Methods */
			std::vector<double> calcCCXDiag(const int &p, const int &q, const CFG *cfg);
			std::vector<double> calcDiag();
			std::vector<double> calcDiag(const WaveFunction *wave_function);
			std::vector<double> calcSigma(const DataFOIS &data_fois,
										  const std::vector<double> &ci_vector,
										  const WaveFunction *wfn_cipsi);

			std::vector<std::vector<double>> calcGuess(const std::vector<double> &diag);

			std::vector<std::unordered_map<std::string, double>>
			mapPreviousCIVector(const std::vector<std::vector<double>> &ci_vectors);

			void createHFConf(std::set<std::string> &cfgs_new,
							  wfn_ptr &wave_function,
							  std::map<int, std::set<std::string>> &sfs,
							  std::map<int, std::map<int, std::string>> &sfs_map__idx_to_sf,
							  std::map<int, std::map<std::string, int>> &sfs_map__sf_to_idx);

			void createAllSFs(const int &nue_max,
							  std::map<int, std::set<std::string>> &sfs,
							  std::map<int, std::map<int, std::string>> &sfs_map__idx_to_sf,
							  std::map<int, std::map<std::string, int>> &sfs_map__sf_to_idx);

			void createAllSFsRecursive(const int &nue, char step,
									   double s, int i, std::string sf,
									   std::map<int, std::set<std::string>> &sfs);

			void readCSFs(const std::string &csfs_fname);

			void runDriver(std::vector<double> &ci_energies,
						   std::vector<std::vector<double>> &ci_vectors,
						   std::vector<std::vector<double>> &energies_per_iter,
						   std::vector<std::unordered_map<std::string, double>> &previous_ci_coeffs_map);

			void solveCI(std::vector<double> &ci_energies,
						 std::vector<std::vector<double>> &ci_vectors,
						 std::vector<std::vector<double>> &energies_per_iter,
						 std::vector<std::unordered_map<std::string, double>> &previous_ci_coeffs_map);

			/* Encapsulated utilities */
			class CIPSI;
			class Connections;
			class CouplingCoeffs;
			class PrefixAlgorithm;
			std::unique_ptr<CIPSI> cipsi;
			std::unique_ptr<Connections> connections;
			std::unique_ptr<CouplingCoeffs> coupling_coeffs;
			std::unique_ptr<PrefixAlgorithm> prefix_algorithm;

			/* Member variables */
			double core_energy;
			double spin;
			int iter_sci;
			int min_nue;
			int multiplicity;
			int n_els;
			int n_orbs;
			int n_roots;

			std::map<int, std::set<std::string>> spin_functions;
			std::map<int, std::map<int, std::string>> sfs_map__idx_to_sf;
			std::map<int, std::map<std::string, int>> sfs_map__sf_to_idx;
			std::map<nonet, cc_map> ccs_2el;
			std::map<quintet, cc_map> ccs_1el;
			std::map<quintet, cc_map> ccs_dia;

			std::set<std::string> cfgs_new;

			vec2d one_el_ham;
			vec2d one_el_ints;
			vec4d two_el_ints;

			std::vector<double> ci_energies;
			std::vector<std::vector<double>> ci_vectors;
			std::vector<std::vector<double>> energies_per_iter;

			std::vector<std::unordered_map<std::string, double>> previous_ci_coeffs_map;

			connection_map_1el connections_1el;
			connection_map_2el connections_2el;
			connection_map_dia connections_dia;

			wfn_ptr wave_function;


			// Temporary timing variables
			double t_solveCI = 0, t_constructCCs = 0, t_connections_construction = 0,
				   t_connections_merge = 0, t_diagonalize = 0;

			double t_selection_setup = 0, t_generateCFGsAndConnections = 0, t_appendSpinFunctions = 0,
				   t_setUpCIPSIWaveFunction = 0, t_constructSFPairs = 0, t_sel_constructCCs = 0,
				   t_doCIPSIPruning = 0, t_appendCFGsAndCSFs = 0, t_selection_total = 0;

			double t_constructCCs1 = 0, t_ccs_1el = 0, t_ccs_2el = 0, t_ccs_dia = 0;
		};
	}
}
