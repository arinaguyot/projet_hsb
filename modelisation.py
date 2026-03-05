import math

def calcul_depot_bloc1_complet(generation_cible, Q_tot_L_min, distribution_gouttes, C0_neb, temps_nebu_min, debit_liquide_mL_min, t_pause_s):
    """
    Calcule la densité de vecteurs sur une bronche (Bloc 1) en respectant 
    le modèle en SÉRIE (les générations précédentes filtrent l'aérosol).
    """
    
    # Dictionnaire de Weibel : {Gen: [n_bronches, D_cm, L_cm, theta_rad]}
    weibel = {
        0: [1, 1.80, 12.0, 0.0],       # Trachée
        1: [2, 1.22, 4.76, 0.52],      # Bronches souches
        2: [4, 0.83, 1.90, 0.52],      # Bronches lobaires
        3: [8, 0.56, 0.76, 0.52],      # Bronches segmentaires
        4: [16, 0.45, 1.27, 0.52],     # Bronches subsegmentaires
        5: [32, 0.35, 1.07, 0.52]      # Petites bronches
    }
    
    if generation_cible not in weibel:
        return "Erreur : Génération non répertoriée."

    # Constantes physiques
    rho_eau = 1000.0  
    mu_air = 1.81e-5  
    g = 9.81          
    Q_tot_m3_s = (Q_tot_L_min * 1e-3) / 60.0

    # ==========================================
    # 1. INITIALISATION À L'ENTRÉE DE LA TRACHÉE
    # ==========================================
    vol_liquide_tot_mL = debit_liquide_mL_min * temps_nebu_min
    vecteurs_initiaux_totaux = vol_liquide_tot_mL * C0_neb
    
    # On crée un dictionnaire pour suivre les vecteurs "survivants" par taille de goutte
    # Format : {diamètre_um: nombre_de_vecteurs_encore_en_vol}
    vecteurs_en_vol = {}
    for d_um, fraction in distribution_gouttes:
        vecteurs_en_vol[d_um] = vecteurs_initiaux_totaux * fraction

    vecteurs_deposes_cible = 0.0 # Ce qu'on veut trouver à la fin

    # ==========================================
    # 2. DESCENTE DANS L'ARBRE (Le modèle en série)
    # ==========================================
    # On boucle de la trachée (0) jusqu'à la génération cible incluse
    for gen in range(generation_cible + 1):
        
        n_bronches, D_cm, L_cm, theta_rad = weibel[gen]
        D = D_cm * 1e-2
        L = L_cm * 1e-2
        
        # Dynamique des fluides locale
        Q_local = Q_tot_m3_s / n_bronches           
        Surface_section = math.pi * (D / 2)**2
        V_air = Q_local / Surface_section           
        t_total_sedimentation = (L / V_air) + t_pause_s

        # On calcule le dépôt pour chaque taille de goutte DANS CETTE GÉNÉRATION
        for d_um in vecteurs_en_vol:
            d_g = d_um * 1e-6
            
            # Calcul des probabilités de dépôt
            V_sed = (rho_eau * d_g**2 * g) / (18 * mu_air) 
            exposant_sed = (4 * V_sed * t_total_sedimentation) / (math.pi * D)
            P_sed = 1 - math.exp(-exposant_sed)
            
            Stk = (rho_eau * d_g**2 * V_air) / (18 * mu_air * D) 
            P_imp = 1 - math.exp(-Stk * theta_rad)
            
            P_tot = P_sed + P_imp - (P_sed * P_imp)
            
            # Nombre de vecteurs de cette taille qui se déposent ICI (sur toute la génération)
            depot_local = vecteurs_en_vol[d_um] * P_tot
            
            # Si c'est notre génération cible, on sauvegarde le résultat !
            if gen == generation_cible:
                vecteurs_deposes_cible += depot_local
            
            # MISE À JOUR CRUCIALE : On retire les vecteurs déposés de ceux qui continuent le voyage
            vecteurs_en_vol[d_um] = vecteurs_en_vol[d_um] - depot_local

    # ==========================================
    # 3. TRADUCTION EN DENSITÉ DE SURFACE (C0) POUR LE BLOC 2
    # ==========================================
    # On a le total déposé sur TOUTES les bronches de la génération cible.
    # On divise par le nombre de bronches pour en avoir une seule.
    vecteurs_sur_UNE_bronche = vecteurs_deposes_cible / weibel[generation_cible][0]
    
    surface_bronche_cm2 = math.pi * weibel[generation_cible][1] * weibel[generation_cible][2]
    densite_surface = vecteurs_sur_UNE_bronche / surface_bronche_cm2
    
    return densite_surface

# --- TEST ---
distribution = [(1.0, 0.20), (3.0, 0.50), (6.0, 0.30)]
C0 = calcul_depot_bloc1_complet(4, 15.0, distribution, 1e10, 30, 0.3, 3.0)
print(f"DENSITÉ INITIALE (C0) corrigée des pertes en amont : {C0:.2e} vg/cm2")