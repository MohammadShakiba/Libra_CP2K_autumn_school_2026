"""fssh_nbra: copy this file, rename it, and edit the selected options.

Keep exactly one active choice in each category. The alternatives below are
those listed in the source tutorial/recipes, not a catalogue of every Libra
parameter. Some require forces, widths, or other data absent in NBRA.

Shared settings follow ../fssh_nbra.py: HST NACs, adiabatic LD integrator 2,
thermal acceptance 31, no forces or momentum rescaling, and no separate state
tracking/phase correction. Method choices come from 9_step4/recipes. Their
mixed 31/32 acceptance and IDA/DISH Hvib integrators are harmonized here so
comparisons use the same propagation and acceptance choices.

Selected: tsh_method=0, decoherence_algo=-1, decoherence_times_type=-1.
DISH uses its own events (no extra damping); revised DISH uses FSSH plus
correction 7. MSDM and both DISH choices need the supplied decoherence rates.
If changing from LD to an Hvib integrator, review tracking and phase correction
together. If changing the NAC algorithm, also rebuild Hvib consistently
(hvib_update_method=1); the reader currently constructs it with HST NACs.
"""

def load(dyn_general):

    # Uncomment one of the options in each of the categories below:
    #====== How to update Hamiltonian ===================
    #dyn_general.update({"ham_update_method":0}) # don't update any Hamiltonians
    #dyn_general.update({"ham_update_method":1})  # recompute only diabatic Hamiltonian, common choice for model Hamiltonians
    dyn_general.update({"ham_update_method":2})  # recompute only adiabatic Hamiltonian; use with file-based or on-the-fly workflows


    #====== How to transform the Hamiltonians between representations ============
    dyn_general.update( {"ham_transform_method":0 }) # don't do any transforms; usually for NBRA or on-the-fly workflows,
                                                      # so you don't override the read values
    #dyn_general.update( {"ham_transform_method":1 }) # diabatic->adiabatic according to internal diagonalization
    #dyn_general.update( {"ham_transform_method":2 }) # diabatic->adiabatic according to internally stored basis transformation matrix
    #dyn_general.update( {"ham_transform_method":3 }) # adiabatic->diabatic according to internally stored basis transformation matrix
    #dyn_general.update( {"ham_transform_method":4 }) # adiabatic->diabatic according to local diabatization method

    #====== How to get the time-overlaps in the dynamics ========
    dyn_general.update( {"time_overlap_method":0 })  # don't update time-overlaps - maybe they are already pre-computed and read
    #dyn_general.update( {"time_overlap_method":1 }) # explicitly compute it from the wavefunction info; common for model Hamiltonians

    #================== How to compute NACs ===============================
    #dyn_general.update({"nac_update_method":0 })  # just read from files
    #dyn_general.update({"nac_update_method":1})  # explicit NAC calculations - let's just focus on this one for now
    dyn_general.update({"nac_update_method":2, "nac_algo":0})  # HST algo
    #dyn_general.update({"nac_update_method":2, "nac_algo":1})  # NPI algo

    #============== How to compute vibronic Hamiltonian ==============
    dyn_general.update( {"hvib_update_method":0 }) # don't update Hvib; maybe because we read it from files
    #dyn_general.update( {"hvib_update_method":1 }) # recompute diabatic and adiabatic Hvib from the Ham and NACs in those reps

    #=========== Ehrenfest or state-resolved options ===========
    # This is what we use with any of the TSH-based methods - in all cases here, we would
    # use "rep_force":1 so that we are guided by the forces derived from the adiabatic surfaces.
    # In Ehrenfest cases though, the forces can be computed using only diabatic properties though
    dyn_general.update( {"force_method":0, "rep_force":1} ) # don't compute forces
    #dyn_general.update( {"force_method":1, "rep_force":1} ) # state-resolved (e.g. TSH) with adiabatic forces
    #dyn_general.update( {"force_method":2, "rep_force":1} ) # for Ehrenfest in adiabatic rep
    #dyn_general.update( {"force_method":2, "rep_force":0} ) # for Ehrenfest in diabatic rep


    #============ Types of surface hopping acceptance and momenta rescaling options =================
    #dyn_general.update({"hop_acceptance_algo":10, "momenta_rescaling_algo":100 })  # accept and rescale based on total energy, do not reverse on frustrated
    #dyn_general.update({"hop_acceptance_algo":10, "momenta_rescaling_algo":101 })  # accept and rescale based on total energy, reverse on frustrated
    #dyn_general.update({"hop_acceptance_algo":20, "momenta_rescaling_algo":200 })  # accept and rescale based on NAC vectors, do not reverse on frustrated
    #dyn_general.update({"hop_acceptance_algo":20, "momenta_rescaling_algo":201 })  # accept and rescale based on NAC vectors, reverse on frustrated
    #dyn_general.update({"hop_acceptance_algo":21, "momenta_rescaling_algo":200 })  # accept and rescale based on force differences, do not reverse on frustrated
    #dyn_general.update({"hop_acceptance_algo":21, "momenta_rescaling_algo":201 })  # accept and rescale based on force differences, reverse on frustrated
    dyn_general.update({"hop_acceptance_algo":31, "momenta_rescaling_algo":0 })  # quantum Boltzmann acceptance; no rescaling
    #dyn_general.update({"hop_acceptance_algo":32, "momenta_rescaling_algo":0 })  # classical Maxwell-Boltzmann acceptance; no rescaling

    #============ Surface hop proposal methodology =================
    #dyn_general.update({"tsh_method":-1 }) # adiabatic dynamics, no hops [ default ]
    dyn_general.update({"tsh_method":0 }) # Fewest Switches Surface Hopping (FSSH)
    #dyn_general.update({"tsh_method":1 }) # Global Flux Surface Hopping (GFSH)
    #dyn_general.update({"tsh_method":2 }) # Markov-State Surface Hopping (MSSH)
    #dyn_general.update({"tsh_method":3, "rep_lz":0 }) # Landau-Zener (LZ) options
    #dyn_general.update({"tsh_method":4, "rep_lz":0 }) # Zhu-Nakamura (ZN) options
    #dyn_general.update({"tsh_method":5 }) # DISH
    #dyn_general.update({"tsh_method":6 }) # MASH
    #dyn_general.update({"tsh_method":7 }) # FSSH2
    #dyn_general.update({"tsh_method":8 }) # FSSH3
    #dyn_general.update({"tsh_method":9 }) # GFSH (original)

    #=========== Decoherence options =================
    dyn_general.update({ "decoherence_algo":-1}) # no (additional) decoherence
    #dyn_general.update({ "decoherence_algo":0}) # msdm
    #dyn_general.update({ "decoherence_algo":1}) # IDA
    #dyn_general.update({ "decoherence_algo":2}) # A-FSSH, not yet ready
    #dyn_general.update({ "decoherence_algo":3}) # BCSH
    #dyn_general.update({ "decoherence_algo":4}) # MF-SD
    #dyn_general.update({ "decoherence_algo":5}) # SHXF
    #dyn_general.update({ "decoherence_algo":6}) # MQCXF
    #dyn_general.update({ "decoherence_algo":7}) # DISH, revised 2023

    #=========== Decoherence times (MSDM, MF-SD, DISH variants) ==================
    # Manual times use the decoherence_rates matrix supplied by the simulation.
    # For Schwartz options, first define A as the inverse-width MATRIX.
    dyn_general.update({"decoherence_times_type":-1 }) # No decoherence times, infinite decoherence times
    #dyn_general.update( { "decoherence_times_type":0 } )  # manual decoherence times
    #dyn_general.update( { "decoherence_times_type":1, "decoherence_C_param": 1.0, "decoherence_eps_param":0.1 } )  # EDC + default params
    #dyn_general.update( { "decoherence_times_type":2, "schwartz_decoherence_inv_alpha":A } ) # Schwartz version 1
    #dyn_general.update( { "decoherence_times_type":3, "schwartz_decoherence_inv_alpha":A } ) # Schwartz version 2

    #======= Various decoherence-related parameters =====================
    dyn_general.update({"dephasing_informed":0})
    #dyn_general.update({"dephasing_informed":1}) # requires decoherence_rates and ave_gaps


    #======= DISH-specific ======================
    #dyn_general.update({"dish_decoherence_event_option":0}) # event at the decoherence time
    dyn_general.update({"dish_decoherence_event_option":1}) # stochastic event sampling

    #=========== Phase correction of SSY =================
    dyn_general.update({"do_ssy":0 }) # no SSY correction
    #dyn_general.update({"do_ssy":1 }) # SSY correction

    #=========== What to integrate ==================
    # solve TD-SE in diabatic representation
    #dyn_general.update({"rep_tdse":0, "electronic_integrator":-1 })   # no propagation
    #dyn_general.update({"rep_tdse":0, "electronic_integrator":0 })    # Lowdin exp_ with 2-point Hvib_dia
    #dyn_general.update({"rep_tdse":0, "electronic_integrator":1 })    # based on QTAG propagator
    #dyn_general.update({"rep_tdse":0, "electronic_integrator":2 })    # based on modified QTAG propagator (Z at two times)
    #dyn_general.update({"rep_tdse":0, "electronic_integrator":3 })    # non-Hermitian integrator with 2-point Hvib_dia

    # solve TD-SE in adiabatic representation
    #dyn_general.update({"rep_tdse":1, "electronic_integrator":-1 })  # no propagation
    #dyn_general.update({"rep_tdse":1, "electronic_integrator":0 })   # ld, with crude splitting,  with exp_
    #dyn_general.update({"rep_tdse":1, "electronic_integrator":1 })   # ld, with symmetric splitting, with exp_
    #dyn_general.update({"rep_tdse":1, "electronic_integrator":2 })   # ld, original, with exp_
    #dyn_general.update({"rep_tdse":1, "electronic_integrator":3 })   # 1-point, Hvib integration, with exp_
    #dyn_general.update({"rep_tdse":1, "electronic_integrator":4 })   # 2-points, Hvib integration, with exp_
    dyn_general.update({"rep_tdse":1, "electronic_integrator":5 })   # 2-points, Hvib, integration with the second-point correction of Hvib, with exp_
    #dyn_general.update({"rep_tdse":1, "electronic_integrator":6 })   # same as 4, but without projection matrices (T_new = I)
    #dyn_general.update({"rep_tdse":1, "electronic_integrator":10 })  # same as 0, but with rotations
    #dyn_general.update({"rep_tdse":1, "electronic_integrator":11 })  # same as 1, but with rotations
    #dyn_general.update({"rep_tdse":1, "electronic_integrator":12 })  # same as 2, but with rotations
    #dyn_general.update({"rep_tdse":1, "electronic_integrator":13 })  # same as 3, but with rotations
    #dyn_general.update({"rep_tdse":1, "electronic_integrator":14 })  # same as 4, but with rotations
    #dyn_general.update({"rep_tdse":1, "electronic_integrator":15 })  # same as 5, but with rotations

    # solve QCLE in diabatic representation
    #dyn_general.update({"rep_tdse":2, "electronic_integrator":0 })  # mid-point Hvib, using exp_

    # solve QCLE in adiabatic representation
    #dyn_general.update({"rep_tdse":3, "electronic_integrator":0 })  # mid-point Ham with the second-point correction of Hvib, using exp_
    #dyn_general.update({"rep_tdse":3, "electronic_integrator":1 })  # using Zhu Liouvillian THIS IS NOT JUST A DIFFERENT INTEGRATOR!!!!
    #dyn_general.update({"rep_tdse":3, "electronic_integrator":10 }) # same as 0 but with rotations

    #=========== Disable state tracking and phase corrections explicitly for the LD integrators ===============
    #State tracking algorithm:
    #  - -1: use LD approach, it includes phase correction too [ default ]
    #  - 0: no state tracking
    #  - 1: method of Kosuke Sato (may fail by getting trapped into an infinite loop)
    #  - 2: Munkres-Kuhn (Hungarian) algorithm
    #  - 21: ChatGPT-generated Munkres-Kuhn (Hungarian) algorithm
    #  - 3: experimental stochastic algorithm, the original version with elimination (known problems)
    #  - 32: experimental stochastic algorithms with all permutations (too expensive)
    #  - 33: the improved stochastic algorithm with good scaling and performance, on par with the mincost
    #  - 4: new, experimental force-based tracking
    #  - 5: SVD-reformulated LD approach of Granucci and Persico
    #dyn_general.update({"state_tracking_algo":-1 })
    dyn_general.update({"state_tracking_algo":0 })
    #dyn_general.update({"state_tracking_algo":1 })
    #dyn_general.update({"state_tracking_algo":2 })
    #dyn_general.update({"state_tracking_algo":21 })
    #dyn_general.update({"state_tracking_algo":3 })
    #dyn_general.update({"state_tracking_algo":32 })
    #dyn_general.update({"state_tracking_algo":33 })
    #dyn_general.update({"state_tracking_algo":4 })
    #dyn_general.update({"state_tracking_algo":5 })


    # Choice of the scaling function for the cost matrix in the Munkres-Kuhn (Hungarian) algorithm
    # Options:
    # - 0 : exp(- alpha^2 * |dE_ij|^2 )  [default]
    # - 1 : exp(- alpha * |dE_ij| )
    # - 2 : exp(-alpha * max(dE_ij, 0) )
    # - anything else:  1  - no scaling
    #dyn_general.update({"MK_scaling_function":0})
    #dyn_general.update({"MK_scaling_function":1})
    #dyn_general.update({"MK_scaling_function":2})
    dyn_general.update({"MK_scaling_function":3}) # no scaling

    # Controlling the range of states involved in state tracking
    # the larger value means narrower range of states is
    # considered when determining possible state tracking effects
    # Doesn't matter for LD integrators
    dyn_general.update({"MK_alpha":1100.0})

    #The algorithm to correct phases on adiabatic states
    #
    #Options:
    #  - 0: no phase correction
    #  - 1: according to our phase correction algorithm [ default ]
    # phase correction doesn't matter, if we use the LD integrator
    dyn_general.update({"do_phase_correction":0 })
    #dyn_general.update({"do_phase_correction":1 })

    #New phase correction, directly applied to NACs. Intended to be used mostly with state_tracking_algo == 4,
    #although can be useful with other state treacking algorithms. Should not be used together with
    #`do_phase_correction`
    #Options:
    #  - 0: no correction [ default ]
    #  - 1: do this correction
    dyn_general.update({"do_nac_phase_correction":0 })
    #dyn_general.update({"do_nac_phase_correction":1 })



