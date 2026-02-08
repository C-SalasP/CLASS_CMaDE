

include("module_CMaDE.jl")
using .CMaDEModule





h_LCDM    = 0.678100 
Ω0k_LCDM  = 0.0
Ω0g_LCDM  = 5.37815e-05
Ω0ur_LCDM = 3.71799e-05
Ω0b_LCDM  = 0.0486773
Ω0cdm_LCDM = 0.261206
Ω0L_LCDM = 1 - Ω0g_LCDM - Ω0ur_LCDM - Ω0b_LCDM - Ω0cdm_LCDM -Ω0k_LCDM


# --- Parámetros cosmológicos para CMaDE---
h        = 0.6934 - 0.0066 + 0.0045                       # Parámetro de Hubble reducido
Ω0k      = -0.0056  + 0.0051                            # Curvatura espacial
Ω0g      = 5.37815e-05                       # Densidad de fotones
Ω0ur     = 3.71799e-05                       # Densidad de radiación ultrarelativista (neutrinos)
Ω0b      =  (0.0223 + 0.0001)/(h*h)                      # Densidad de materia bariónica
Ω0DM     =  (0.2984 + 0.0060 ) - Ω0b      #0.271206                           # Densidad de materia oscura fría (CDM)
Ω0de = 1.0 - Ω0g - Ω0ur - Ω0b - Ω0DM - Ω0k     # Asumiendo universo plano (Ωk ≈ 0)

if (Ω0de < 0)
    error("La densidad de energía oscura calculada es negativa. Por favor, revise los parámetros cosmológicos proporcionados.")
end

Q_coupling= -0.063 + 0.0309

k_c = 1.0

has_CMaDE= "yes"
save_perturbations = "no"
background = "no"
gauge = "Newtonian"
matter_source_in_current_gauge_CMaDE = "yes"
matter_source_in_current_gauge_LCDM = "yes"




z_pk=0.0

status=CMaDEModule.run_CMaDE(Ω0g, Ω0ur, Ω0b, Ω0DM, Ω0de, h, k_c, Ω0k, z_pk, Q_coupling, has_CMaDE,save_perturbations, background, matter_source_in_current_gauge_CMaDE, gauge)

if status == 0

CMaDEModule.run_LCDM(Ω0g_LCDM, Ω0ur_LCDM, Ω0b_LCDM, Ω0cdm_LCDM, Ω0L_LCDM, h_LCDM, Ω0k_LCDM, z_pk, save_perturbations, background, matter_source_in_current_gauge_LCDM, gauge)

CMaDEModule.plot_latest_pk(z_pk)
CMaDEModule.plot_latest_cl()
#CMaDEModule.plot_background()
CMaDEModule.run_CMaDE_Gauge(Ω0g, Ω0ur, Ω0b, Ω0DM, Ω0de, h, k_c, Ω0k, z_pk, Q_coupling, has_CMaDE,save_perturbations, background, matter_source_in_current_gauge_CMaDE)

CMaDEModule.plot_gauges_CMaDE(z_pk)
end


