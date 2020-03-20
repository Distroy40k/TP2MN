echo "Ceci est le programme de test de gemm"
echo "V(X) = vecteur avec pour chaque coordonées la valeur X"
echo "M(X) = Matrice remplie de X"
echo "Le calcul est décrit avant de lancer le test, sous la forme"
echo "alpha * A * B + beta * C; taille des matrices ; mode"
echo "###############################"
echo "###############################"
read -p "Press enter to continue"

##############################
## TEST s mode ##
###############################
echo
echo
echo "Calcul: V(1) * M(1) * M(0) + V(1) * M(1); 2x2 ; s"
echo "Res attendu :  M(1)"
echo "######################"
read -p "Press enter to continue"
./test_gemm s fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/m_2_2_0.txt fichiers/Rbeta_1.txt fichiers/m_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(0) + V(2) * M(1); 2x2 ; s"
echo "Res attendu :  M(2)"
echo "######################"
read -p "Press enter to continue"
./test_gemm s fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/m_2_2_0.txt fichiers/Rbeta_2.txt fichiers/m_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(1) + V(1) * M(1); 2x2 ; s"
echo "Res attendu :  M(3)"
echo "######################"
read -p "Press enter to continue"
./test_gemm s fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/m_2_2_1.txt fichiers/Rbeta_1.txt fichiers/m_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(1) + V(0) * M(1); 2x2 ; s"
echo "Res attendu :  M(2)"
echo "######################"
read -p "Press enter to continue"
./test_gemm s fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/m_2_2_1.txt fichiers/Rbeta_0.txt fichiers/m_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(0) + V(0) * M(1); 2x2 ; s"
echo "Res attendu :  M(0)"
echo "######################"
read -p "Press enter to continue"
./test_gemm s fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/m_2_2_0.txt fichiers/Rbeta_0.txt fichiers/m_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(1) + V(0) * M(1); 5x5 ; s"
echo "Res attendu :  M(5)"
echo "######################"
read -p "Press enter to continue"
./test_gemm s fichiers/Ralpha_1.txt fichiers/m_5_5_1.txt fichiers/m_5_5_1.txt fichiers/Rbeta_0.txt fichiers/m_5_5_1.txt 5

##############################
## TEST D mode ##
###############################
echo
echo
echo "Calcul: V(1) * M(1) * M(0) + V(1) * M(1); 2x2 ; s"
echo "Res attendu :  M(1)"
echo "######################"
read -p "Press enter to continue"
./test_gemm d fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/m_2_2_0.txt fichiers/Rbeta_1.txt fichiers/m_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(0) + V(2) * M(1); 2x2 ; s"
echo "Res attendu :  M(2)"
echo "######################"
read -p "Press enter to continue"
./test_gemm d fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/m_2_2_0.txt fichiers/Rbeta_2.txt fichiers/m_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(1) + V(1) * M(1); 2x2 ; s"
echo "Res attendu :  M(3)"
echo "######################"
read -p "Press enter to continue"
./test_gemm d fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/m_2_2_1.txt fichiers/Rbeta_1.txt fichiers/m_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(1) + V(0) * M(1); 2x2 ; s"
echo "Res attendu :  M(2)"
echo "######################"
read -p "Press enter to continue"
./test_gemm d fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/m_2_2_1.txt fichiers/Rbeta_0.txt fichiers/m_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(0) + V(0) * M(1); 2x2 ; s"
echo "Res attendu :  M(0)"
echo "######################"
read -p "Press enter to continue"
./test_gemm d fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/m_2_2_0.txt fichiers/Rbeta_0.txt fichiers/m_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(1) + V(0) * M(1); 5x5 ; s"
echo "Res attendu :  M(5)"
echo "######################"
read -p "Press enter to continue"
./test_gemm d fichiers/Ralpha_1.txt fichiers/m_5_5_1.txt fichiers/m_5_5_1.txt fichiers/Rbeta_0.txt fichiers/m_5_5_1.txt 5



##############################
## TEST C mode ##
###############################
echo
echo
echo "Calcul: V(1) * M(1) * M(0) + V(1) * M(1); 2x2; c"
echo "Res attendu :  M(1)"
echo "######################"
read -p "Press enter to continue"
./test_gemm c fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/mC_2_2_0.txt fichiers/Cbeta_1.txt fichiers/mC_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(1) + V(1) * M(1); 2x2; c"
echo "Res attendu :  Matrice de coefficients égaux à 1 suivi de 5 (partie réelle puis imaginaire)"
echo "######################"
read -p "Press enter to continue"
./test_gemm c fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/mC_2_2_1.txt fichiers/Cbeta_1.txt fichiers/mC_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(1) + V(0) * M(1); 2x2; c"
echo "Res attendu :  Matrice de coefficients égaux à 0 suivi de 4 (partie réelle puis imaginaire)"
echo "######################"
read -p "Press enter to continue"
./test_gemm c fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/mC_2_2_1.txt fichiers/Cbeta_0.txt fichiers/mC_2_2_1.txt 2
echo
echo
echo "Calcul: V(2) * M(1) * M(1) + V(0) * M(1); 2x2; c"
echo "Res attendu :  Matrice de coefficients égaux à 0 suivi de 8 (partie réelle puis imaginaire)"
echo "######################"
read -p "Press enter to continue"
./test_gemm c fichiers/Calpha_2.txt fichiers/mC_2_2_1.txt fichiers/mC_2_2_1.txt fichiers/Cbeta_0.txt fichiers/mC_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(0) + V(0) * V(1); 2x2; c"
echo "Res attendu :  M(0)"
echo "######################"
read -p "Press enter to continue"
./test_gemm c fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/mC_2_2_0.txt fichiers/Cbeta_0.txt fichiers/mC_2_2_1.txt 2
echo
echo
echo "Calcul: V(2) * M(1) * M(1) + V(0) * M(1); 5x5; c"
echo "Res attendu :  Matrice de coefficients égaux à 0 suivi de 20 (partie réelle puis imaginaire)"
echo "######################"
read -p "Press enter to continue"
./test_gemm c fichiers/Calpha_2.txt fichiers/mC_5_5_1.txt fichiers/mC_5_5_1.txt fichiers/Cbeta_0.txt fichiers/mC_5_5_1.txt 5



##############################
## TEST Z mode ##
###############################
echo
echo
echo "Calcul: V(1) * M(1) * M(0) + V(1) * M(1); 2x2; c"
echo "Res attendu :  M(1)"
echo "######################"
read -p "Press enter to continue"
./test_gemm z fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/mC_2_2_0.txt fichiers/Cbeta_1.txt fichiers/mC_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(1) + V(1) * M(1); 2x2; c"
echo "Res attendu :  Matrice de coefficients égaux à 1 suivi de 5 (partie réelle puis imaginaire)"
echo "######################"
read -p "Press enter to continue"
./test_gemm z fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/mC_2_2_1.txt fichiers/Cbeta_1.txt fichiers/mC_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(1) + V(0) * M(1); 2x2; c"
echo "Res attendu :  Matrice de coefficients égaux à 0 suivi de 4 (partie réelle puis imaginaire)"
echo "######################"
read -p "Press enter to continue"
./test_gemm z fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/mC_2_2_1.txt fichiers/Cbeta_0.txt fichiers/mC_2_2_1.txt 2
echo
echo
echo "Calcul: V(2) * M(1) * M(1) + V(0) * M(1); 2x2; c"
echo "Res attendu :  Matrice de coefficients égaux à 0 suivi de 8 (partie réelle puis imaginaire)"
echo "######################"
read -p "Press enter to continue"
./test_gemm z fichiers/Calpha_2.txt fichiers/mC_2_2_1.txt fichiers/mC_2_2_1.txt fichiers/Cbeta_0.txt fichiers/mC_2_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * M(0) + V(0) * V(1); 2x2; c"
echo "Res attendu :  M(0)"
echo "######################"
read -p "Press enter to continue"
./test_gemm z fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/mC_2_2_0.txt fichiers/Cbeta_0.txt fichiers/mC_2_2_1.txt 2
echo
echo
echo "Calcul: V(2) * M(1) * M(1) + V(0) * M(1); 5x5; c"
echo "Res attendu :  Matrice de coefficients égaux à 0 suivi de 20 (partie réelle puis imaginaire)"
echo "######################"
read -p "Press enter to continue"
./test_gemm z fichiers/Calpha_2.txt fichiers/mC_5_5_1.txt fichiers/mC_5_5_1.txt fichiers/Cbeta_0.txt fichiers/mC_5_5_1.txt 5