#### MODELE AEROSOL

#---------Bloc 1----------

import math

def calcul_depot_bloc1(generation, Q_tot_L_min, d_goutte_um, C0_neb, temps_nebu_min, debit_liquide_mL_min):
    """
    Calcule la densité de vecteurs déposés sur le mucus d'une bronche.
    Basé sur le modèle de Rudolf et al. (sédimentation + impaction).
    
    Paramètres :
    - generation : Numéro de la génération bronchique (ex: 4 pour bronches subsegmentaires)
    - Q_tot_L_min : Débit d'air total du respirateur EVLP (L/min)
    - d_goutte_um : Diamètre de la gouttelette d'aérosol (micromètres)
    - C0_neb : Concentration de vecteurs dans le nébuliseur (vg/mL)
    - temps_nebu_min : Durée totale de la nébulisation (minutes)
    - debit_liquide_mL_min : Débit de liquide consommé par le nébuliseur (mL/min)
    """
    
    ### Anatomie (modèle de Wiebel)
    # Format : {Génération: [Nombre_bronches, Diamètre_cm, Longueur_cm, Angle_bifurcation_rad]}
    weibel = {
        0: [1, 1.80, 12.0, 0.0],       # Trachée (Tube EVLP)
        1: [2, 1.22, 4.76, 0.52],      # Bronches souches (Angle ~30 deg)
        2: [4, 0.83, 1.90, 0.52],      # Bronches lobaires
        3: [8, 0.56, 0.76, 0.52],      # Bronches segmentaires
        4: [16, 0.45, 1.27, 0.52],     # Bronches subsegmentaires
        5: [32, 0.35, 1.07, 0.52]      # Petites bronches
    }
    
    if generation not in weibel:
        return "Erreur : Génération non répertoriée dans ce dictionnaire."
        
    n_bronches, D_cm, L_cm, theta_rad = weibel[generation]
    
    ### Conversions 
    d_g = d_goutte_um * 1e-6
    D = D_cm * 1e-2
    L = L_cm * 1e-2
    Q_tot_m3_s = (Q_tot_L_min * 1e-3) / 60.0
    
    ### Constantes physiques 
    rho_eau = 1000.0  # kg/m3 (densité de la goutte)
    mu_air = 1.81e-5  # Pa.s (viscosité de l'air)
    g = 9.81          # m/s2 (gravité)

    # ==========================================
    # 2. DYNAMIQUE DES FLUIDES DANS LA BRONCHE
    # ==========================================
    Q_local = Q_tot_m3_s / n_bronches           # Débit dans UNE seule bronche
    Surface_section = math.pi * (D / 2)**2
    V_air = Q_local / Surface_section           # Vitesse de l'air (m/s)

    ###Probabilités de déposition (modèle de Rudolf)
    
    # A. Sédimentation (Chute par gravité)
    V_sed = (rho_eau * d_g**2 * g) / (18 * mu_air) # Vitesse de chute terminale
    exposant_sed = (4 * V_sed * L) / (math.pi * D * V_air)
    P_sed = 1 - math.exp(-exposant_sed)
    
    # B. Impaction (Inertie dans les virages)
    Stk = (rho_eau * d_g**2 * V_air) / (18 * mu_air * D) # Nombre de Stokes
    P_imp = 1 - math.exp(-Stk * theta_rad)
    
    # C. Probabilité totale (Sédimentation OU Impaction)
    P_tot = P_sed + P_imp - (P_sed * P_imp)

    # ==========================================
    # 4. TRADUCTION EN THÉRAPIE GÉNIQUE (Lien Bloc 2)
    # ==========================================
    # Volume total de liquide entrant dans CETTE bronche spécifique
    vol_liquide_tot_mL = debit_liquide_mL_min * temps_nebu_min
    vol_liquide_local_mL = vol_liquide_tot_mL / n_bronches
    
    # Nombre de vecteurs qui entrent puis qui se déposent
    vecteurs_entrants = vol_liquide_local_mL * C0_neb
    vecteurs_deposes = vecteurs_entrants * P_tot
    
    # Calcul de la densité de surface (vg/cm2)
    surface_bronche_cm2 = math.pi * D_cm * L_cm
    densite_surface = vecteurs_deposes / surface_bronche_cm2
    
    return densite_surface, P_sed, P_imp, P_tot


# ==========================================
# EXEMPLE D'UTILISATION (TEST DU CODE)
# ==========================================

# Paramètres de test (à adapter à votre protocole EVLP)
generation_cible = 4          # Bronches subsegmentaires
debit_EVLP = 15.0             # L/min
taille_goutte = 3.0           # Micromètres (typique d'un nébuliseur standard)
titre_viral = 1e10            # vg/mL dans la cuve du nébuliseur
temps_traitement = 30         # minutes
debit_nebuliseur = 0.3        # mL/min

# Exécution de la fonction
densite, prob_sed, prob_imp, prob_tot = calcul_depot_bloc1(
    generation_cible, debit_EVLP, taille_goutte, titre_viral, temps_traitement, debit_nebuliseur
)

# Affichage des résultats
print(f"--- RÉSULTATS BLOC 1 (Génération {generation_cible}) ---")
print(f"Probabilité Sédimentation : {prob_sed*100:.2f} %")
print(f"Probabilité Impaction     : {prob_imp*100:.2f} %")
print(f"Probabilité Totale        : {prob_tot*100:.2f} %")
print("-" * 35)
print(f"DENSITÉ INITIALE SUR LE MUCUS (C0) : {densite:.2e} vg/cm2")