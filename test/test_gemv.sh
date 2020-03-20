echo "Ceci est le programme de test de gemv"
echo "V(X) = vecteur avec pour chaque coordonées la valeur X"
echo "M(X) = Matrice remplie de X"
echo "Le calcul est décrit avant de lancer le test, sous la forme"
echo "alpha * A * X + beta * Y; taille_du_vecteur; mode"
echo "###############################"
echo "###############################"
read -p "Press enter to continue"

##############################
## TEST s mode ##
###############################
echo
echo
echo "Calcul: V(1) * M(1) * V(0) + V(1) * V(1); 2; s"
echo "Res attendu :  v(1)"
echo "######################"
read -p "Press enter to continue"
./test_gemv s fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/x_2_0.txt fichiers/Rbeta_1.txt fichiers/y_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * V(0) + V(2) * V(1); 2; s"
echo "Res attendu :  v(2)"
echo "######################"
read -p "Press enter to continue"
./test_gemv s fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/x_2_0.txt fichiers/Rbeta_2.txt fichiers/y_2_1.txt 2
echo "Calcul: V(1) * M(1) * V(1) + V(1) * V(1); 2; s"
echo "Res attendu :  v(3)"
echo "######################"
read -p "Press enter to continue"
./test_gemv s fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/x_2_1.txt fichiers/Rbeta_1.txt fichiers/y_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * V(1) + V(0) * V(1); 2; s"
echo "Res attendu :  v(2)"
echo "######################"
read -p "Press enter to continue"
./test_gemv s fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/x_2_1.txt fichiers/Rbeta_0.txt fichiers/y_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * V(0) + V(0) * V(1); 2; s"
echo "Res attendu :  v(0)"
echo "######################"
read -p "Press enter to continue"
./test_gemv s fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/x_2_0.txt fichiers/Rbeta_0.txt fichiers/y_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * V(1) + V(0) * V(1); 128; s"
echo "Res attendu :  v(128)"
echo "######################"
read -p "Press enter to continue"
./test_gemv s fichiers/Ralpha_1.txt fichiers/m_128_128_1.txt fichiers/x_128_1.txt fichiers/Rbeta_0.txt fichiers/y_128_1.txt 128

##############################
## TEST D mode ##
###############################
echo
echo
echo "Calcul: V(1) * M(1) * V(0) + V(1) * V(1); 2; d"
echo "Res attendu :  v(1)"
echo "######################"
read -p "Press enter to continue"
./test_gemv d fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/x_2_0.txt fichiers/Rbeta_1.txt fichiers/y_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * V(0) + V(2) * V(1); 2; d"
echo "Res attendu :  v(2)"
echo "######################"
read -p "Press enter to continue"
./test_gemv d fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/x_2_0.txt fichiers/Rbeta_2.txt fichiers/y_2_1.txt 2
echo
echo "Calcul: V(1) * M(1) * V(1) + V(1) * V(1); 2; d"
echo "Res attendu :  v(3)"
echo "######################"
read -p "Press enter to continue"
./test_gemv d fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/x_2_1.txt fichiers/Rbeta_1.txt fichiers/y_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * V(1) + V(0) * V(1); 2; d"
echo "Res attendu :  v(2)"
echo "######################"
read -p "Press enter to continue"
./test_gemv d fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/x_2_1.txt fichiers/Rbeta_0.txt fichiers/y_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * V(0) + V(0) * V(1); 2; d"
echo "Res attendu :  v(0)"
echo "######################"
read -p "Press enter to continue"
./test_gemv d fichiers/Ralpha_1.txt fichiers/m_2_2_1.txt fichiers/x_2_0.txt fichiers/Rbeta_0.txt fichiers/y_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * V(1) + V(0) * V(1); 128; d"
echo "Res attendu :  v(128)"
echo "######################"
read -p "Press enter to continue"
./test_gemv d fichiers/Ralpha_1.txt fichiers/m_128_128_1.txt fichiers/x_128_1.txt fichiers/Rbeta_0.txt fichiers/y_128_1.txt 128

##############################
## TEST C mode ##
###############################
echo
echo
echo "Calcul: V(1) * M(1) * V(0) + V(1) * V(1); 2; c"
echo "Res attendu :  0 2"
echo "               0 2"
echo "######################"
read -p "Press enter to continue"
./test_gemv c fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/xC_2_0.txt fichiers/Cbeta_1.txt fichiers/yC_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * V(1) + V(1) * V(1); 2; c"
echo "Res attendu :  -4 6"
echo "               -4 6"
echo "######################"
read -p "Press enter to continue"
./test_gemv c fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/xC_2_1.txt fichiers/Cbeta_1.txt fichiers/yC_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * V(1) + V(0) * V(1); 2; c"
echo "Res attendu :  -4 4"
echo "               -4 4"
echo "######################"
read -p "Press enter to continue"
./test_gemv c fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/xC_2_1.txt fichiers/Cbeta_0.txt fichiers/yC_2_1.txt 2
echo
echo
echo "Calcul: V(2) * M(1) * V(1) + V(0) * V(1); 2; c"
echo "Res attendu :  -8 8"
echo "               -8 8"
echo "######################"
read -p "Press enter to continue"
./test_gemv c fichiers/Calpha_2.txt fichiers/mC_2_2_1.txt fichiers/xC_2_1.txt fichiers/Cbeta_0.txt fichiers/yC_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * V(0) + V(0) * V(1); 2; c"
echo "Res attendu :  0 0"
echo "               0 0"
echo "######################"
read -p "Press enter to continue"
./test_gemv c fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/xC_2_0.txt fichiers/Cbeta_0.txt fichiers/yC_2_1.txt 2
echo
echo
echo "Calcul: V(2) * M(1) * V(1) + V(0) * V(1); 128; c"
echo "Res attendu :  -512 512"
echo "               512 512"
echo "               ... ..."
echo "######################"
read -p "Press enter to continue"
./test_gemv c fichiers/Calpha_2.txt fichiers/mC_128_128_1.txt fichiers/xC_128_1.txt fichiers/Cbeta_0.txt fichiers/yC_128_1.txt 128



##############################
## TEST Z mode ##
###############################
echo
echo
echo "Calcul: V(1) * M(1) * V(0) + V(1) * V(1); 2; z"
echo "Res attendu :  0 2"
echo "               0 2"
echo "######################"
read -p "Press enter to continue"
./test_gemv z fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/xC_2_0.txt fichiers/Cbeta_1.txt fichiers/yC_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * V(1) + V(1) * V(1); 2; z"
echo "Res attendu :  -4 6"
echo "               -4 6"
echo "######################"
read -p "Press enter to continue"
./test_gemv z fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/xC_2_1.txt fichiers/Cbeta_1.txt fichiers/yC_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * V(1) + V(0) * V(1); 2; z"
echo "Res attendu :  -4 4"
echo "               -4 4"
echo "######################"
read -p "Press enter to continue"
./test_gemv z fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/xC_2_1.txt fichiers/Cbeta_0.txt fichiers/yC_2_1.txt 2
echo
echo
echo "Calcul: V(2) * M(1) * V(1) + V(0) * V(1); 2; z"
echo "Res attendu :  -8 8"
echo "               -8 8"
echo "######################"
read -p "Press enter to continue"
./test_gemv z fichiers/Calpha_2.txt fichiers/mC_2_2_1.txt fichiers/xC_2_1.txt fichiers/Cbeta_0.txt fichiers/yC_2_1.txt 2
echo
echo
echo "Calcul: V(1) * M(1) * V(0) + V(0) * V(1); 2; z"
echo "Res attendu :  0 0"
echo "               0 0"
echo "######################"
read -p "Press enter to continue"
./test_gemv z fichiers/Calpha_1.txt fichiers/mC_2_2_1.txt fichiers/xC_2_0.txt fichiers/Cbeta_0.txt fichiers/yC_2_1.txt 2
echo
echo
echo "Calcul: V(2) * M(1) * V(1) + V(0) * V(1); 128; z"
echo "Res attendu :  -512 512"
echo "               512 512"
echo "               ... ..."
echo "######################"
read -p "Press enter to continue"
./test_gemv z fichiers/Calpha_2.txt fichiers/mC_128_128_1.txt fichiers/xC_128_1.txt fichiers/Cbeta_0.txt fichiers/yC_128_1.txt 128
