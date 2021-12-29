import spd_analyze


dir_name_agafia = '!SPDStandResults//AGAFIA//21_06_08__18_54'
dir_name_bronislav = '!SPDStandResults//BRONISLAV//21_07_20__21_03'
dir_name_alyonushka = "!SPDStandResults//Alyonushka//21_12_02__13_47"

spd_agafia = spd_analyze.SPD(dir_name_agafia)
spd_agafia = spd_agafia.get_spd()


spd_bronislav = spd_analyze.SPD(dir_name_bronislav)
spd_bronislav = spd_bronislav.get_spd()

spd_alyonuushka = spd_analyze.SPD(dir_name_alyonushka)
spd_alyonuushka = spd_alyonuushka.get_spd()

