declare -a List_pat=(
AAAS_13071969
AAMV_27011994
AGBB_26121972
ASCS_31121976
BBMA_19041977
BCUV_02091982
C_DG_15111995
CACB_27091965
CAMO_08111988
CARM_14011974
CDPR_25031986
CEBT_30091989
CJGV_24021957
CMCO_17051975
CMLL_11051985
COGC_18022000
CRCM_04051980
DFCM_04071987
EERV_01041955
EFHB_26121994
ELFN_09021961
FABN_19111984
FFGV_27111965
FJAM_10012001
FRCA_12081956
GADM_08111983
GAQJ_21061987
GJBA_03051990
IABF_16091997
IABM_03061982
IASA_27071974
J_BF_06111999
JAPV_03111995
JAZV_27081991
JCZD_03071963
JDRP_06081956
JFFA_04101976
JFSO_05061967
JGRF_24091991
JMRR_10111976
L_GM_04061960
LAAB_29101982
LDEV_31011988
LMRG_21041962
MAEM_17051957
MCAG_07091991
MFCM_21011984
MFSM_02021973
MGRN_22101979
MJCS_29091978
MLAD_09021957
MLCB_11041981
MPMM_22111977
MTVV_27111974
NLJA_20011972
OAPC_01041962
RABG_09091981
RALF_14041986
RAPA_19121977
RASG_07061964
REOO_10011989
RFCR_12111978
RSVT_17091955
SAAA_10051991
SAJC_22021972
SAVS_23061979
UENM_22121979
VAPH_11101985
YVSH_02051983#68
     )
	
	

		 from=0
		 to=${#List_pat[@]}
		 
		 arguments=($@)
		 index=0
		 na=${#arguments[@]}
		 while [ ${index} -lt ${na} ]; do
		 	AR=${arguments[index]}
		 	case ${AR} in
	 			--to=*)
	 			to=${AR/*=/""}
	 			;;
 			    --from=*)
 			    from=${AR/*=/""}
 			    ;;				
		 	esac
	
		 	index=$((index + 1))
		 done	
	



		 DB=DBNC_03/neuroCOVID
		 NUBE="/Volumes/GoogleDrive-103235447575506129142/Mi unidad" 
		 
		 
		 
		 
		 for i in `seq ${from} ${to}`
		 do

		     for CHECK in EEG EYE
			 do

			 if [ ! -d "/Volumes/DBNC_03/neuroCOVID/DATA/${List_pat[$i]}/${CHECK}" ] && [ -d "/Volumes/DBNC_03/neuroCOVID/DATA/${List_pat[$i]}/${CHECK}" ]; 			then	 
				echo '/Volumes/DBNC_03/neuroCOVID/DATA/'${List_pat[$i]}'/'${CHECK}
				cp -R -f  "/Volumes/DBNC_03/neuroCOVID/DATA/${List_pat[$i]}/${CHECK}"  "/Volumes/DBNC_03/neuroCOVID/DATA/${List_pat[$i]}/${CHECK}" 
		     fi
		
		 done
		 done