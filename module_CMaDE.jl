# ==========================================================
# module_CMaDE.jl
# Módulo con funciones del modelo CMaDE
# ==========================================================
module CMaDEModule

export run_CMaDE
export run_LCDM
export plot_latest_pk
export plot_latest_cl
export plot_gauges_CMaDE
export plot_background
export run_CMaDE_Gauge


using DelimitedFiles
using PyPlot
using Glob
using Statistics
using Interpolations
using DifferentialEquations
using Printf
using DelimitedFiles
using PyPlot
using Glob
using Statistics
using PyCall   

# Importa pyplot
plt = pyimport("matplotlib.pyplot")

# Importa GridSpec desde matplotlib
gridspec = pyimport("matplotlib.gridspec")


# ----------------------------------------------------------
# Sistema de ecuaciones diferenciales
# ----------------------------------------------------------
function omega_system_loga!(dO, Os, p, loga)
    rho_g, rho_ur, rho_b, rho_dm, rho_de = Os
    Ω0k, k_c, Q_coupling = p
    a = exp(loga)
    
    
    H_squared = rho_g + rho_ur + rho_b + rho_dm + rho_de + Ω0k / a^2
    if(H_squared < 1e-8)
       error("H^2 es negativo o muy pequeño, lo que puede causar errores numéricos.")
    end
    H = sqrt(H_squared)
    drho_de = Q_coupling*sqrt(6) * rho_de^(3/2) / (a * H * π)

    dO[1] = -4 * rho_g
    dO[2] = -4 * rho_ur
    dO[3] = -3 * rho_b
    dO[4] = -3 * rho_dm - k_c*drho_de
    dO[5] = drho_de


end

# ----------------------------------------------------------
# Función principal que resuelve, crea el .ini y ejecuta CLASS
# ----------------------------------------------------------
function run_CMaDE(Ω0g, Ω0ur, Ω0b, Ω0DM, Ω0de, h, k_c, Ω0k, z_pk, Q_coupling, has_CMaDE::String,save_perturbations::String, save_background::String, matter_source_in_current_gauge::String, gauge::String)

    Ω0 = [Ω0g, Ω0ur, Ω0b, Ω0DM, Ω0de]
    p = (Ω0k, k_c, Q_coupling)
    aspan = (0.0, log(1e-14))
    status = 0
    sol = nothing
    try
        prob = ODEProblem(omega_system_loga!, Ω0, aspan, p)
        sol = solve(prob, Rodas5(), reltol=1e-18, abstol=1e-20)  
    catch e
        println(" Error al resolver las ecuaciones diferenciales.")
        println(" Detalles del error: ", e)
        return 1
    end

    rho_g  = sol[1, :]
    rho_ur = sol[2, :]
    rho_b  = sol[3, :]
    rho_dm = sol[4, :]
    rho_de = sol[5, :]
    rho_ini = (rho_g[end], rho_ur[end], rho_b[end], rho_dm[end], rho_de[end])

   #println("densidad de energía inicial CMaDE de: ", rho_ini[5])
   #println("densidad de energía inicial CMaDE dm: ", rho_ini[4])

    ini_text = """
    has_CMaDE = $has_CMaDE
    save_perturbations = $save_perturbations
    modes = s
    ic = ad

    k_c = $(k_c)
    Q_coupling = $(Q_coupling)
    # --- Densidades iniciales ---
    Omega_ini_CMaDE_dm = $(rho_ini[4])
    Omega_ini_CMaDE_de = $(rho_ini[5])
  


    # --- Densidades actuales ---
    Omega_g  = $(Ω0g)
    Omega_ur = $(Ω0ur)
    Omega_b  = $(Ω0b)
    Omega_k= $(Ω0k)

    matter_source_in_current_gauge = $matter_source_in_current_gauge
    # --- Parámetros adicionales ---
    h = $(h)
    output = tCl ,mPk
    z_pk = $(z_pk)
    background_verbose = 0
    gauge = $gauge
    write_background = $save_background
    """

    open("AutomaticCMaDE.ini", "w") do f
        write(f, ini_text)
    end

    try
        run(`./class AutomaticCMaDE.ini`)
        println("se ejecutó CLASS para CMaDE.")
    catch e
        println(" Error al ejecutar CLASS para CMaDE: ")
        println(" Detalles del error: ", e)
        status = 1
    end
    return status
end 





function run_LCDM(Ω0g, Ω0ur, Ω0b, Ω0DM, Ω0de, h, Ω0k, z_pk, save_perturbations::String, save_background::String, matter_source_in_current_gauge::String, gauge::String)
    ini_text2 = """
    save_perturbations = $save_perturbations
    modes = s
    ic = ad


    # --- Densidades actuales ---
    Omega0_g  = $(Ω0g)
    Omega0_ur = $(Ω0ur)
    Omega0_b  = $(Ω0b)
    Omega_k= $(Ω0k)
    Omega_L = $(Ω0de)
    Omega_cdm = $(Ω0DM)


    matter_source_in_current_gauge = $matter_source_in_current_gauge

    # --- Parámetros adicionales ---
    h = $(h)
    output = tCl,mPk
    z_pk = $(z_pk)
    background_verbose = 0
    gauge = $(gauge)
    write_background = $save_background
    """


    open("AutomaticLCDM.ini", "w") do f
        write(f, ini_text2)
    end

    try
        run(`./class AutomaticLCDM.ini`)
                println("se ejecutó CLASS para LCDM.")
    catch e
        println(" Error al ejecutar CLASS para LCDM: ", e)
    end

end 








function plot_latest_pk(z_plot)
    lcdm_files = sort(glob("output/AutomaticLCDM*_pk.dat"))
    cmade_files = sort(glob("output/AutomaticCMaDE*_pk.dat"))

    if isempty(lcdm_files) || isempty(cmade_files)
        println("No se encontraron archivos de pk para LCDM o CMaDE")
        return
    end

    lcdm_file = lcdm_files[end]
    cmade_file = cmade_files[end]

    lcdm_data = readdlm(lcdm_file, skipstart=5)
    cmade_data = readdlm(cmade_file, skipstart=5)

    # --- Definir variables dentro de la función ---
    k_lcdm = lcdm_data[:,1]
    P_lcdm = lcdm_data[:,2]

    k_cmade = cmade_data[:,1]
    P_cmade = cmade_data[:,2]

    # --- Recortar a la tercera parte en adelante ---
    N_lcdm = length(k_lcdm)
    N_cmade = length(k_cmade)

    idx_lcdm = 40 : N_lcdm
    idx_cmade = 40 : N_cmade


    k_lcdm = k_lcdm[idx_lcdm]
    P_lcdm = P_lcdm[idx_lcdm]

    k_cmade = k_cmade[idx_cmade]
    P_cmade = P_cmade[idx_cmade]


    # --- Interpolación para comparar ---
    # Interpolamos CMaDE en los mismos k que LCDM
    itp = LinearInterpolation(k_cmade, P_cmade, extrapolation_bc=Line())
    P_cmade_interp = itp.(k_lcdm)

    # --- Calcular error absoluto ---
    error_abs = abs.(P_lcdm .- P_cmade_interp)./ abs.(P_cmade_interp)

    # --- Graficar ---


    fig = figure(figsize=(10,6))
    # Creamos dos subplots con height_ratios diferentes
    # gs = GridSpec(nrows, ncols, height_ratios=[...])
    gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])

    ax1 = subplot(gs[1, 0])          # subplot superior
    ax2 = subplot(gs[2, 0], sharex=ax1)  # subplot inferior comparte eje x

    # Subplot 1
    ax1.plot(k_lcdm, P_lcdm, label="LCDM", lw=2, linestyle="--")    
    ax1.plot(k_cmade, P_cmade, label="CMaDE")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_ylabel("P(k) [(Mpc)^3]")
    ax1.set_title("Matter Power Spectrum at z=$z_plot")
    ax1.legend()
    ax1.grid(false)

    # Subplot 2
    ax2.plot(k_lcdm, error_abs, color="red", lw=2)
    ax2.set_xlabel("k [1/Mpc]")
    ax2.set_ylabel("Error relativo")
    #ax2.set_title("Absolute error")
    ax2.grid(false)

    tight_layout()

    savefig("zz_MisGraficasCLASS/Pk.png",dpi=500)
    show()


end

function plot_latest_cl()
    # --- Buscar los archivos más recientes ---
    lcdm_files = sort(glob("output/AutomaticLCDM*_cl.dat"))  # todos los archivos LCDM
    cmade_files = sort(glob("output/AutomaticCMaDE*_cl.dat")) # todos los archivos CMaDE

    if isempty(lcdm_files) || isempty(cmade_files)
        println("No se encontraron archivos de Cl para LCDM o CMaDE")
        return
    end

    lcdm_file = lcdm_files[end]  # el más reciente
    cmade_file = cmade_files[end]

    # --- Leer datos ---
    lcdm_data = readdlm(lcdm_file, skipstart=12)   # ajustar skipstart según tu archivo
    cmade_data = readdlm(cmade_file, skipstart=12)

    # Asumiendo columnas: 1 = l, 2 = C_l^TT
    l_lcdm = lcdm_data[:,1]
    Cl_lcdm = lcdm_data[:,2]

    l_cmade = cmade_data[:,1]
    Cl_cmade = cmade_data[:,2]

    # --- Interpolación para comparar ---
    itp = LinearInterpolation(l_cmade, Cl_cmade, extrapolation_bc=Line())
    Cl_cmade_interp = itp.(l_lcdm)

    # --- Calcular error relativo ---
    error_rel = abs.(Cl_lcdm .- Cl_cmade_interp) ./ abs.(Cl_cmade_interp)

    # --- Graficar ---
    fig = figure(figsize=(10,6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])

    ax1 = subplot(gs[1, 0])
    ax2 = subplot(gs[2, 0], sharex=ax1)

    # Subplot 1: C_l TT
    ax1.plot(l_lcdm, Cl_lcdm, label="LCDM", lw=2)
    ax1.plot(l_cmade, Cl_cmade, label="CMaDE", lw=2, linestyle="--")

    #ax1.set_xscale("log")
    #ax1.set_yscale("log")
    ax1.set_ylabel("C_l^TT")
    ax1.set_title("CMB TT Power Spectrum")
    ax1.legend()
    ax1.grid(false)

    # Subplot 2: error relativo
    ax2.plot(l_lcdm, error_rel, color="red", lw=2)
    ax2.set_xlabel("Multipole l")
    ax2.set_ylabel("Error relativo")
    ax2.grid(false)

    tight_layout()
    savefig("zz_MisGraficasCLASS/CMB.png",dpi=500)
    show()
end

function plot_background()
      # --- Buscar los archivos más recientes ---
    lcdm_files = sort(glob("output/AutomaticLCDM*_background.dat"))  # todos los archivos LCDM
    cmade_files = sort(glob("output/AutomaticCMaDE*_background.dat")) # todos los archivos CMaDE

    if isempty(lcdm_files) || isempty(cmade_files)
        println(" No se encontraron archivos de pk para LCDM o CMaDE")
        return
    end

    lcdm_file = lcdm_files[end]  
    cmade_file = cmade_files[end]

    # --- Leer datos ---
    lcdm_data = readdlm(lcdm_file, skipstart=5)  
    cmade_data = readdlm(cmade_file, skipstart=5)

    a_cmade = cmade_data[:, 1] 
    Ω_r_CMaDE = (cmade_data[:, 10] .+ cmade_data[:, 15]) ./ (cmade_data[:, 5].^2) 
    Ω_dm_CMaDE =cmade_data[:, 12] ./ (cmade_data[:, 5].^2)
    Ω_de_CMaDE = cmade_data[:, 13] ./ (cmade_data[:, 5].^2) 
    Ω_b_CMaDE = cmade_data[:, 11] ./ (cmade_data[:, 5].^2) 

    a_lcdm = lcdm_data[:, 1] 
    Ω_r_LCDM = (lcdm_data[:, 10] .+ lcdm_data[:, 14]) ./ lcdm_data[:, 16] 
    Ω_dm_LCDM =lcdm_data[:, 12] ./ lcdm_data[:, 16] 
    Ω_de_LCDM = lcdm_data[:, 13] ./ lcdm_data[:, 16] 
    Ω_b_LCDM = lcdm_data[:, 11] ./ lcdm_data[:, 16]
    


    figure(figsize=(10, 7))

    plot(a_lcdm, Ω_dm_LCDM, label="Ω_dm (LCDM)", color="cyan", linestyle="-", linewidth=3)
    plot(a_cmade, Ω_dm_CMaDE, label="Ω_dm (CMaDE)", color="blue", linestyle="--", linewidth=3)

    plot(a_lcdm, Ω_b_LCDM, label="Ω_b (LCDM)", color="green", linestyle="-", linewidth=3)
    plot(a_cmade, Ω_b_CMaDE, label="Ω_b (CMaDE)", color="lime", linestyle="--", linewidth=3)

    plot(a_lcdm, Ω_de_LCDM, label="Ω_de (LCDM)", color="red", linestyle="-", linewidth=3)
    plot(a_cmade, Ω_de_CMaDE, label="Ω_de (CMaDE)", color="black", linestyle="--", linewidth=3)

    plot(a_lcdm, Ω_r_LCDM, label="Ω_r (LCDM)", color="orange", linestyle="-", linewidth=3)
    plot(a_cmade, Ω_r_CMaDE, label="Ω_r (CMaDE)", color="yellow", linestyle="--", linewidth=3)

    title("Comparación de Evolución de Densidades")
    xscale("log")
    xlabel(" (a)")
    ylabel("Ω(a)")
    legend(loc="best")
    grid(false)
    show()
end











function plot_gauges_CMaDE(z_plot)
    sync_files  = sort(glob("output/AutomaticGaugeCMaDE*_pk.dat"))
    cmade_files = sort(glob("output/AutomaticCMaDE*_pk.dat"))

    if isempty(sync_files) || isempty(cmade_files)
        println("No se encontraron archivos de pk para gauge síncrono o CMaDE")
        return
    end

    sync_file  = sync_files[end]
    cmade_file = cmade_files[end]

    sync_data  = readdlm(sync_file,  skipstart=5)
    cmade_data = readdlm(cmade_file, skipstart=5)

    # --- Definir variables dentro de la función ---
    k_sync = sync_data[:,1]
    P_sync = sync_data[:,2]

    k_cmade = cmade_data[:,1]
    P_cmade = cmade_data[:,2]

    # --- Recortar a la tercera parte en adelante ---
    N_sync  = length(k_sync)
    N_cmade = length(k_cmade)

    idx_sync  = 1 : N_sync
    idx_cmade = 1 : N_cmade

    k_sync = k_sync[idx_sync]
    P_sync = P_sync[idx_sync]

    k_cmade = k_cmade[idx_cmade]
    P_cmade = P_cmade[idx_cmade]

    # --- Interpolación para comparar ---
    # Interpolamos CMaDE en los mismos k que gauge síncrono
    itp = LinearInterpolation(k_cmade, P_cmade, extrapolation_bc=Line())
    P_cmade_interp = itp.(k_sync)

    # --- Calcular error relativo ---
    error_abs = abs.(P_sync .- P_cmade_interp) ./ abs.(P_cmade_interp)

    # --- Graficar ---
    fig = figure(figsize=(10,6))

    gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])

    ax1 = subplot(gs[1, 0])
    ax2 = subplot(gs[2, 0], sharex=ax1)

    # Subplot 1
    ax1.plot(k_sync,  P_sync,  label="Synchronous gauge", lw=2, linestyle="--")
    ax1.plot(k_cmade, P_cmade, label="Newtonian gauge", lw=2)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_ylabel("P(k) [(Mpc)^3]")
    ax1.set_title("Matter Power Spectrum at z=$z_plot")
    ax1.legend()
    ax1.grid(false)

    # Subplot 2
    ax2.plot(k_sync, error_abs, lw=2)
    ax2.set_xlabel("k [1/Mpc]")
    ax2.set_ylabel("Error relativo")
    ax2.grid(false)

    tight_layout()
    savefig("zz_MisGraficasCLASS/Pk_sync_vs_newtonian.png", dpi=500)
    show()
end










function run_CMaDE_Gauge(Ω0g, Ω0ur, Ω0b, Ω0DM, Ω0de, h, k_c, Ω0k, z_pk, Q_coupling, has_CMaDE::String,save_perturbations::String, save_background::String, matter_source_in_current_gauge::String)

    Ω0 = [Ω0g, Ω0ur, Ω0b, Ω0DM, Ω0de]
    p = (Ω0k, k_c, Q_coupling)
    aspan = (0.0, log(1e-14))
    status = 0
    sol = nothing
    try
        prob = ODEProblem(omega_system_loga!, Ω0, aspan, p)
        sol = solve(prob, Rodas5(), reltol=1e-18, abstol=1e-20)  
    catch e
        println(" Error al resolver las ecuaciones diferenciales.")
        println(" Detalles del error: ", e)
        return 1
    end

    rho_g  = sol[1, :]
    rho_ur = sol[2, :]
    rho_b  = sol[3, :]
    rho_dm = sol[4, :]
    rho_de = sol[5, :]
    rho_ini = (rho_g[end], rho_ur[end], rho_b[end], rho_dm[end], rho_de[end])

   


    ini_text3 = """
    has_CMaDE = $has_CMaDE
    save_perturbations = $save_perturbations
    modes = s
    ic = ad

    k_c = $(k_c)
    Q_coupling = $(Q_coupling)
    # --- Densidades iniciales ---
    Omega_ini_CMaDE_dm = $(rho_ini[4])
    Omega_ini_CMaDE_de = $(rho_ini[5])


    # --- Densidades actuales ---
    Omega_g  = $(Ω0g)
    Omega_ur = $(Ω0ur)
    Omega_b  = $(Ω0b)
    Omega_k= $(Ω0k)

    matter_source_in_current_gauge = $(matter_source_in_current_gauge)
    # --- Parámetros adicionales ---
    h = $(h)
    output = tCl ,mPk
    z_pk = $(z_pk)
    background_verbose = 0
    gauge = synchronous
    write_background = $save_background
    """

    open("AutomaticGaugeCMaDE.ini", "w") do f
        write(f, ini_text3)
    end

    try
        run(`./class AutomaticGaugeCMaDE.ini`)
        println("se ejecutó CLASS para CMaDE en el gauge synchronous.")
    catch e
        println(" Error al ejecutar CLASS para CMaDE en el gauge synchronous: ")
        println(" Detalles del error: ", e)
        status = 1
    end
    return status
end 





end # module