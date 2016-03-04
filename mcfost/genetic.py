def mcfost_genetic(grid_name, parameter_file, fitting_routine,
    simu_config_file=None,root_dir=None, n_nodes=20,
    omp_num_threads=8,walltime=24.0,
    waiting_time=30.,options=None,n_models=100,
    pmut=0.2,pmut_max=0.5,pmut_min=1e-2,generation_max=100,elitism=1,
    no_mcfost_compute=False,progressive_plot=True, rt=True)  :
    """mcfost genetic algorithm main routine



    """
