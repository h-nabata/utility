#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N 260 // the maximum number of charactors in a line (the number of bytes) !!!! ## if N is more than 256, this program does NOT work @ IMS (9501)
#define L 400 // the maximum number of charactor array  ## reduced @ IMS (1628)
#define K 400 // the maximum number of charactor array  ## reduced @ IMS (1628)
#define MAX 5000 // the maximum number of atoms in the system  ## reduced @ IMS (1628)  ## increased @ IMS (3208)

/// ### EQ composition and structure sorting tool ### ///

/*
my memo

 gcc EQCS_ver1.5.c -o EQCS_ver1.5.out -lm -O2

* Use JOBNAME constantly
* 

# Options
> crystal mode is automatically aplied.  ;for TV system
> option: "scaled output" "-o"   ;EQリストをスケーリングする
> option: "composition sort" "-c"   ;EQとパスを組成別に分類し、LUP計算しやすいようにファイリングする
** 吐き出すファイルは"EQCS_[JOBNAME]_xxx.log"とする

> atomname_TVexcluded には[1]から原子名を格納している（[0]の中身は空）

> 変数 babel_compositionnum は廃止し、inchi表記に統一する（cansmiでは水素結合の配向まで区別してしまうため）
> Iseq はinchi_composition_groupnum に格納している（babelによるsortはこれに従っている）

> LUP,TS・PTに対応(9425 planned)
> LUPの端点、TS・PTの端点を分類したい

> DC_list に対応したい(9511 planned)

#9715 CIFファイルを生成できるように改良
#9821 DOS計算用のinpを生成するように改良
#1628 確保するメモリの量を削減（IMS用に調整） --> 機能してない
#2320 "the same Iseq connection" についても障壁の高さを表示
#3606 原子数の数え方を改善
#4316 PTについて最安定EQを基準とするバリアの値も出力するように変更

*/

////// original functions　// "strcount(char *str1, char *str2)" search *str2 in *str1  // equal to strcmp
double radius(char *);
int digit(int);
int strcount(char *, char *);
int mystrcmp(char *, char *);
void lntrim(char *);
int Zvalue(char *);
double Atomicmass(char *);
//////

int main(int argc, char* argv[]) {

//////***** values about this software *****//////
  char software_version[N] = "1.5c" ; char developed_date[N] = "2024/03/16";
//////*****                            *****//////

/// variables for file processing
  FILE *fp,*output,*scaledEQlist,*rearrangedEQlist,*supercellxyz,*spEQfile,*gjf,*cif,*ciflist,*supercellgjf,*supercellEQlist,*cellVDlist,*temp_babel,*energy_sort,*energy_sort_xyz,*density_sort,*forDOSinp,*forDOSsh,*cytoscapeEQlist,*cytoscapeTSlist,*cytoscapePTlist; // file structure for EQlist sorting
  FILE *sublog,*for_test,*re_output,*re_sublog,*temp_inputfile,*temp_outputfile; // file structure for EQcomposition sorting
  char FILENAME_argument[N], original_EQlist_fname[N], JOBNAME[N], xyzfilename[N], tempfname[N], scaled_EQlist_filename[N], rearranged_EQlist_filename[N], cell_VDlist_filename[N], scaled_cell_VDlist_filename[N], cell_dlist_filename[N], dlogfilename[N], elogfilename[N], exyzfilename[N], spfilename[N], gjffilename[N], ciffilename[N], ciflistfilename[N], forDOSinpfilename[N], forDOSshfilename[N], cytofilename[N]; // file name
  char filename[N], sublogname[N], reinput_filename[N], reoutput_filename[N];
  char str[N], temp_atomname1[MAX][N], temp_atomname2[MAX][N], temp_atomname_for_rearrange[MAX][N], atomname_TVexcluded[MAX][N], temp_buff[N], temp_buff2[N], command[L];
  double unit_cell_mass = 0.00, temp_dnum1, temp_dnum2, temp_dnum3;
  double TSTS_MAXenergy = -9999999.9, PTTS_MAXenergy = -9999999.9;
  double TVextension = 20.0; // default value
  ///###/// identify JOBtype automatically and put a flag
  int JOBtype = 0;
  //# 1:FREQ 2:MIN 3:SADDLE 4:IRC  8:RESTRUCT   13:LUP 16:SC-AFIR 17:DS-AFIR
  ///###///
  int crystal_surface_mode = 0, composition_sorting_option = 0, scaled_EQlist_output_option = 0, rearranged_EQlist_output_option = 0, supercellgen_option = 0, gjfgen_option = 0, cifgen_option = 0, cifbohr_option = 0, forDOSinpgen_option = 0, cytoscapelistgen_option = 0, TVextension_option = 0, temp_option_num = 1, initial_FILENAME = 0; // !! temp_option_num begins from "2" !!
  int no_babel_mode = 0;
  int one_by_one_mode = 0;
  int energy_sort_mode = 0;
  int density_sort_mode = 0;
  int logfile_mode = 0;
  int EQlist2_mode = 0; // if EQ_list_2.log is opened, set "EQlist2_mode = 1"
  int TSlist2_mode = 0; // if TS_list_2.log is opened, set "TSlist2_mode = 1"
  int LUP_mode = 0; // if artificial force is NOT applied , set "LUP_mode = 1"
  ///###///
  int atoms = 0, meta_atoms=0, TVcounting = 0, TVnum = 0, totalTVnum, TVatom_linenum[2], original_atom_order[MAX], rearranged_atom_order[MAX], totalEQinfonum = 0; // count the number of atoms
  int totalelementnamenum = 0, addnewatom = 1, each_atomnum, elementnamelength, elementnamelength1, elementnamelength2; // moved 9821
  int forDOSinp_Species_num_list[MAX], forDOSinp_Species_order_list[MAX], DOSk1 = 4, DOSk2 = 4, DOSk3 = 4, DOScutoff = 200, DOSnp = 8; // added 9821
  int atomcounting = 0, tempEQnum = 0, totalEQnum = 0;
  int tempnum = 1, tempnum2 = 0, slength;
  int framescale_x = 0, framescale_y = 0, framescale_z = 0, xyzframescale_x = 0, xyzframescale_y = 0, xyzframescale_z = 0, sp_scale_factor = 0;
///

/// variables for EQ list up
  int i = 0, j, k, l, m, n, jj, kk, ll, mm;
  int EQinfonum = 0, nmodenum = 0, temp_scaled_atomnum = 0, temp_nmodenum = 0, SPtemp_nmodenum = 0, scaled_nmodenum = 0, SPscaled_nmodenum = 0, totalatoms_digitnum, totalEQ_digitnum;
///

//////*****  initial check for arguments (options: c n e d f l o r s h sp gjf cif cs)
  char asterisk[N] = "***********************************************************************\n";
  char WQ[3] = "\"";
  char usage[5*N] = "\n Options\n   : -[c/C] composition sorting with Open Babel\n   : -[n/N] not using Open Babel (with -[c/C] option)\n   : -[e/E] sort energy\n   : -[d/D] sort density\n   : -[arr/r] generate rearranged EQlist\n   : -[sp/SP *] generate super cells in the direction of the shortest TV\n   : -[gjf] generate gjf files for each EQ\n   : -[cif] generate cif files (in ang)\n   : -[cifb] generate cif files in bohr\n   : -[dos/DOS] generate .inp file\n   : -[dosk/DOSk a b c] specify the k-grid values\n   : -[dosm/DOSm *] specify the mesh cutoff value\n   : -[x/X/tvx/TVx *] extend the TV length along the a axis\n   : -[cs/cyto] generate files for Cytoscape\n   : -[f/F] read the list_2.log file\n   : -[l/L] read the .log file\n   : -[o/O] convert EQ xyz files one by one\n   : -[s/S i j k] generate i*j*k scaled EQ list\n (d, o, s, sp, and dos options only work for the system including TV)\n";
  if( argc < 2 || argc > 18 ) {
    printf( "%s Usage : \"%s [JOBNAME] [option]\" (%s)%s%s", asterisk, argv[0], software_version, usage, asterisk );
    return -1;
  }
  if( strcmp( argv[1] , "-help" ) == 0 || strcmp( argv[1] , "-Help" ) == 0 || strcmp( argv[1] , "--help" ) == 0 || strcmp( argv[1] , "--Help" ) == 0 || strcmp( argv[1] , "--h" ) == 0 || strcmp( argv[1] , "--H" ) == 0 || strcmp( argv[1] , "-h" ) == 0 || strcmp( argv[1] , "-H" ) == 0 ){
    printf( "%s Usage : \"%s [JOBNAME] [option]\" (%s)%s%s", asterisk, argv[0], software_version, usage, asterisk );
    return -1;
  }else{
    sprintf(FILENAME_argument, "%s", argv[1]); // there is room for reconsideration for option input system...
  }
  if( argc > 2 ) { // options
    while ( temp_option_num < argc ){
      if( strcmp( argv[temp_option_num] , "-c" ) == 0 || strcmp( argv[temp_option_num] , "-C" ) == 0 ){
        composition_sorting_option = 1;
      }else if( strcmp( argv[temp_option_num] , "-o" ) == 0 || strcmp( argv[temp_option_num] , "-O" ) == 0 ){
        one_by_one_mode = 1;
      }else if( strcmp( argv[temp_option_num] , "-n" ) == 0 || strcmp( argv[temp_option_num] , "-N" ) == 0 ){
        no_babel_mode = 1;
      }else if( strcmp( argv[temp_option_num] , "-e" ) == 0 || strcmp( argv[temp_option_num] , "-E" ) == 0 ){
        energy_sort_mode = 1;
      }else if( strcmp( argv[temp_option_num] , "-d" ) == 0 || strcmp( argv[temp_option_num] , "-D" ) == 0 ){
        density_sort_mode = 1;
      }else if( strcmp( argv[temp_option_num] , "-f" ) == 0 || strcmp( argv[temp_option_num] , "-F" ) == 0 ){
        EQlist2_mode = 1;
      }else if( strcmp( argv[temp_option_num] , "-l" ) == 0 || strcmp( argv[temp_option_num] , "-L" ) == 0 ){
        logfile_mode = 1;
      }else if( strcmp( argv[temp_option_num] , "-help" ) == 0 || strcmp( argv[temp_option_num] , "-Help" ) == 0 || strcmp( argv[temp_option_num] , "--help" ) == 0 || strcmp( argv[temp_option_num] , "--Help" ) == 0 || strcmp( argv[temp_option_num] , "--h" ) == 0 || strcmp( argv[temp_option_num] , "--H" ) == 0 || strcmp( argv[temp_option_num] , "-h" ) == 0 || strcmp( argv[temp_option_num] , "-H" ) == 0 ){
        printf( "%s Usage : \"%s [JOBNAME] [option]\" (%s)%s%s", asterisk, argv[0], software_version, usage, asterisk );
        return -1;
      }else if( strcmp( argv[temp_option_num] , "-rearr" ) == 0 || strcmp( argv[temp_option_num] , "-arr" ) == 0 || strcmp( argv[temp_option_num] , "-Arr" ) == 0  || strcmp( argv[temp_option_num] , "-R" ) == 0 || strcmp( argv[temp_option_num] , "-r" ) == 0){
        rearranged_EQlist_output_option = 1;
      }else if( strcmp( argv[temp_option_num] , "-s" ) == 0 || strcmp( argv[temp_option_num] , "-S" ) == 0 ){
        scaled_EQlist_output_option = 1;
        if( argc < temp_option_num + 4 ){
          printf( "%s Usage : \"%s [JOBNAME] [option]\" (%s)%s%s", asterisk, argv[0], software_version, usage, asterisk );
          return -1;
        }else{
          temp_option_num++;
          xyzframescale_x = atoi(argv[temp_option_num]);
          temp_option_num++;
          xyzframescale_y = atoi(argv[temp_option_num]);
          temp_option_num++;
          xyzframescale_z = atoi(argv[temp_option_num]);
          framescale_x = xyzframescale_x, framescale_y = xyzframescale_y, framescale_z = xyzframescale_z;
        }
      }else if( strcmp( argv[temp_option_num] , "-sp" ) == 0 || strcmp( argv[temp_option_num] , "-Sp" ) == 0 || strcmp( argv[temp_option_num] , "-SP" ) == 0 ){ // added 9613
        if( argc < temp_option_num + 1 ){
          sp_scale_factor = 1;
        }else{
          temp_option_num++;
          sp_scale_factor = atoi(argv[temp_option_num]);
          if( sp_scale_factor < 1 ){
            sp_scale_factor = 1;
          }else if( sp_scale_factor > 10 ){
            printf("Do you set the super cell scaling factor = %d? ( yes / no )\n >", sp_scale_factor );
            scanf("%s", str);
            if( strcmp( str , "yes" ) == 0 || strcmp( str , "y" ) == 0 || strcmp( str , "Y" ) == 0 ){
              
            }else{
              sp_scale_factor = 1;
            }
          }
        }
        supercellgen_option = 1;
      }else if( strcmp( argv[temp_option_num] , "-gjf" ) == 0 || strcmp( argv[temp_option_num] , "-Gjf" ) == 0 || strcmp( argv[temp_option_num] , "-GJF" ) == 0 ){ // added 9613
        gjfgen_option = 1;
      }else if( strcmp( argv[temp_option_num] , "-cif" ) == 0 || strcmp( argv[temp_option_num] , "-Cif" ) == 0 || strcmp( argv[temp_option_num] , "-CIF" ) == 0 || strcmp( argv[temp_option_num] , "-cifb" ) == 0 || strcmp( argv[temp_option_num] , "-Cifb" ) == 0 || strcmp( argv[temp_option_num] , "-CIFB" ) == 0 ){ // added 9613
        cifgen_option = 1;
        if( strcmp( argv[temp_option_num] , "-cifb" ) == 0 || strcmp( argv[temp_option_num] , "-Cifb" ) == 0 || strcmp( argv[temp_option_num] , "-CIFB" ) == 0){
          cifbohr_option = 1;
        };
      }else if( strcmp( argv[temp_option_num] , "-dos" ) == 0 || strcmp( argv[temp_option_num] , "-Dos" ) == 0 || strcmp( argv[temp_option_num] , "-DOS" ) == 0 ){ // added 9613
        forDOSinpgen_option = 1;
        // そのうち.inpファイルを読み込んで参照できるように改良するかも(9821考え中)
      }else if( strcmp( argv[temp_option_num] , "-dosk" ) == 0 || strcmp( argv[temp_option_num] , "-Dosk" ) == 0 || strcmp( argv[temp_option_num] , "-DOSk" ) == 0 ){ // added 9613
        if( argc < temp_option_num + 4 ){
          printf( "%s Usage : \"%s [JOBNAME] [option]\" (%s)%s%s", asterisk, argv[0], software_version, usage, asterisk );
          return -1;
        }else{
          temp_option_num++;
          DOSk1 = atoi(argv[temp_option_num]);
          temp_option_num++;
          DOSk2 = atoi(argv[temp_option_num]);
          temp_option_num++;
          DOSk3 = atoi(argv[temp_option_num]);
        }
      }else if( strcmp( argv[temp_option_num] , "-dosm" ) == 0 || strcmp( argv[temp_option_num] , "-Dosm" ) == 0 || strcmp( argv[temp_option_num] , "-DOSm" ) == 0 || strcmp( argv[temp_option_num] , "-cutoff" ) == 0 || strcmp( argv[temp_option_num] , "-dosc" ) == 0 ){ // added 9613
        if( argc < temp_option_num + 1 ){
          printf( "%s Usage : \"%s [JOBNAME] [option]\" (%s)%s%s", asterisk, argv[0], software_version, usage, asterisk );
          return -1;
        }else{
          temp_option_num++;
          DOScutoff = atoi(argv[temp_option_num]);
        }
      }else if( strcmp( argv[temp_option_num] , "-dosnp" ) == 0 || strcmp( argv[temp_option_num] , "-Dosnp" ) == 0 || strcmp( argv[temp_option_num] , "-DOSnp" ) == 0 || strcmp( argv[temp_option_num] , "-np" ) == 0 || strcmp( argv[temp_option_num] , "-dosp" ) == 0 ){ // added 9613
        if( argc < temp_option_num + 1 ){
          printf( "%s Usage : \"%s [JOBNAME] [option]\" (%s)%s%s", asterisk, argv[0], software_version, usage, asterisk );
          return -1;
        }else{
          temp_option_num++;
          DOSnp = atoi(argv[temp_option_num]);
        }
      }else if( strcmp( argv[temp_option_num] , "-cs" ) == 0 || strcmp( argv[temp_option_num] , "-cyto" ) == 0 || strcmp( argv[temp_option_num] , "-Cyto" ) == 0 || strcmp( argv[temp_option_num] , "-CYTO" ) == 0 ){ // added 9613
        cytoscapelistgen_option = 1;
      }else if( strcmp( argv[temp_option_num] , "-x" ) == 0 || strcmp( argv[temp_option_num] , "-X" ) == 0 || strcmp( argv[temp_option_num] , "-tvx" ) == 0 || strcmp( argv[temp_option_num] , "-TVX" ) == 0 ){ // added 9613
        TVextension_option = 1;
        if( argc < temp_option_num + 1 ){
          printf( "The default value of TVextension is 20.0 [ang]\n" );
        }else if( strstr( argv[temp_option_num+1] , "-" ) != NULL ){
          temp_option_num++;
          TVextension = atof(argv[temp_option_num]);
        }
      }/*else{
        printf( "  Unknown option \"%s\"\n  Make sure of your options...\n", argv[temp_option_num] );
      }*/
      temp_option_num++;
    }
  }

// determine JOBNAME
  if( strstr( FILENAME_argument , "_EQ_list.log" ) != NULL ){
    slength = strlen(FILENAME_argument);
    sprintf(JOBNAME, "%.*s", slength - 12, FILENAME_argument);
  }else if( strstr( FILENAME_argument , "_EQ_list_2.log" ) != NULL ){
    EQlist2_mode = 1;
    slength = strlen(FILENAME_argument);
    sprintf(JOBNAME, "%.*s", slength - 14, FILENAME_argument);
  }else if( strstr( FILENAME_argument , ".com" ) != NULL || strstr( FILENAME_argument , ".log" ) != NULL ){
    slength = strlen(FILENAME_argument);
    sprintf(JOBNAME, "%.*s", slength - 4, FILENAME_argument);
  }else{
    sprintf(JOBNAME,"%s",FILENAME_argument);
  }

// determine the original file name
  if(EQlist2_mode == 0 && logfile_mode == 0) {
    sprintf( original_EQlist_fname,"%s_EQ_list.log",JOBNAME);
    fp = fopen( original_EQlist_fname,"r"); // open the original file to read
    if(fp == NULL) { // failed to open, return NULL
      sprintf( original_EQlist_fname,"%s_EQ_list_2.log",JOBNAME);
      fp = fopen( original_EQlist_fname,"r"); // open the original file to read
      if(fp == NULL) { // failed to open, return NULL
        printf("%s_EQ_list.log and list_2.log doesn't exist!\n", JOBNAME);
        return -1;
      }else{
        EQlist2_mode = 1;
      }
    }
  }else if(EQlist2_mode == 1 && logfile_mode == 0) {
    sprintf( original_EQlist_fname,"%s_EQ_list_2.log",JOBNAME);
    fp = fopen( original_EQlist_fname,"r"); // open the original file to read
    if(fp == NULL) { // failed to open, return NULL
      printf("%s_EQ_list_2.log doesn't exist!\n", JOBNAME);
      return -1;
    }
  }else if( logfile_mode == 1) { // does not read EQ_list, but .log file
    sprintf( original_EQlist_fname,"%s.log",JOBNAME);
    fp = fopen( original_EQlist_fname,"r"); // open the original file to read
    if(fp == NULL) { // failed to open, return NULL
      printf("%s.log doesn't exist!\n", JOBNAME);
      return -1;
    }
  }

//////*****  initial check for composition
    while( fgets(str, N, fp) != NULL ){
      // atomcounting = 1: counting atom num phase
      // atomcounting = 2: counting EQ num phase
      // atomcounting = 3: finish counting
      if( atomcounting < 2 ){
        if( (strstr(str, "# Geometry of " ) != NULL || strstr(str, "# ITR. " ) != NULL || strstr(str, "#　NODE " ) != NULL || strstr(str, "INITIAL STRUCTURE" ) != NULL) && atomcounting == 0 ){
          atomcounting = 1;
        }
        if( strstr(str, "#" ) == NULL && (strstr(str, "=" ) != NULL || strstr(str, "Energy" ) != NULL || strstr(str, "Threshold" ) != NULL) && atomcounting == 1 ){
          atomcounting = 2;
          TVcounting = 1;
          atoms = atoms - 1; // !! "atoms" doesn't included TVs, "meta_atoms" is the number of atoms including TVs !! //
          meta_atoms = meta_atoms - 1;
          totalEQinfonum++;
        }
        if( strstr(str, "# STEP " ) != NULL && JOBtype == 17 && atomcounting == 1 ){
          atomcounting = 2;
          TVcounting = 1;
          atoms = atoms - 3; // !! "atoms" doesn't included TVs, "meta_atoms" is the number of atoms including TVs !! //
          meta_atoms = meta_atoms - 3;
          totalEQinfonum++;
        }
        // TVでない場合
        if( (atomcounting == 1 && strstr(str, "TV" ) == NULL) && (atomcounting == 1 && strstr(str, "Tv" ) == NULL) ){
          if( strstr(str, "#" ) == NULL ){
            sscanf(str, "%s    %lf    %lf    %lf\n",temp_atomname1[atoms], &temp_dnum1, &temp_dnum2, &temp_dnum3);
            sprintf(temp_atomname2[atoms],"%sX",temp_atomname1[atoms]);
            forDOSinp_Species_num_list[atoms] = Zvalue(temp_atomname2[atoms]);
          }
          atoms++; meta_atoms++;
        }
        // TVの場合
        if( (TVcounting == 0 && strstr(str, "TV" ) != NULL) || (TVcounting == 0 && strstr(str, "Tv" ) != NULL) ){
          crystal_surface_mode = 1;
          TVatom_linenum[TVnum] = meta_atoms-1;
          TVnum++; meta_atoms++;
        }
      }
      if( atomcounting == 2 && strstr(str, ":" ) == NULL ){ // count EQ info line
        totalEQinfonum++;
      }
      if( atomcounting == 2 && strstr(str, ":" ) != NULL ){
        atomcounting = 3;
      }
      if( strstr(str, "# Geometry of EQ" ) != NULL ){
        totalEQnum++;
      }
      if( strstr(str, "FREQ:" ) != NULL && atomcounting == 0 ){ // logfile_mode
        JOBtype = 1;
      }
      if( strstr(str, "MIN:" ) != NULL && atomcounting == 0 ){ // logfile_mode
        JOBtype = 2;
      }
      if( strstr(str, "SADDLE:" ) != NULL && atomcounting == 0 ){ // logfile_mode
        JOBtype = 3;
      }
      if( strstr(str, "IRC:" ) != NULL && atomcounting == 0 ){ // logfile_mode
        JOBtype = 4;
      }
      if( strstr(str, "RESTRUCT:" ) != NULL && atomcounting == 0 ){
        JOBtype = 8;
      }
      if( strstr(str, "LUP:" ) != NULL && atomcounting == 0 ){ // logfile_mode (or normal)
        JOBtype = 13;
      }
      if( strstr(str, "SC-AFIR:" ) != NULL && atomcounting == 0 ){
        JOBtype = 16;
      }
      if( strstr(str, "DS-AFIR:" ) != NULL && atomcounting == 0 ){
        JOBtype = 17;
      }
    }
    totalatoms_digitnum = digit(atoms);
    totalEQ_digitnum = digit(totalEQnum);
    
    for(k = 0;k < meta_atoms;k++){
      unit_cell_mass = unit_cell_mass + Atomicmass(temp_atomname2[k]);
    }
    i = 0;
    
    if( rearranged_EQlist_output_option == 1 ){
      for(m = 0;m < atoms+1 ;m++){
        rearranged_atom_order[m] = m;
        sprintf( temp_atomname_for_rearrange[m], "%s", temp_atomname2[m]);
      }
      for(k = 1;k < atoms+1 ;k++){
        for(m = k+1;m < atoms+1 ;m++){
          if( Zvalue(temp_atomname_for_rearrange[k]) > Zvalue(temp_atomname_for_rearrange[m]) ){
            i = rearranged_atom_order[k];
            sprintf( temp_buff, "%s", temp_atomname_for_rearrange[k]);
            rearranged_atom_order[k] = rearranged_atom_order[m];
            sprintf( temp_atomname_for_rearrange[k], "%s", temp_atomname_for_rearrange[m]);
            rearranged_atom_order[m] = i;
            sprintf( temp_atomname_for_rearrange[m], "%s", temp_buff);
            
          }
        }
      }
      printf("  rearranged atom order\n");
      for(k = 1;k < atoms+1 ;k++){
        for(m = 1;m < atoms+1 ;m++){
          if( k == rearranged_atom_order[m] ){
            original_atom_order[k] = m;
            break;
          }
        }
      }
      for(k = 1;k < atoms+1 ;k++){
        printf("    %-*d :\t%-*d %s \t%d\n", totalatoms_digitnum, k, totalatoms_digitnum, rearranged_atom_order[k], temp_atomname_for_rearrange[k], original_atom_order[k]);
      }
//      return;
    }


    i = 0;
    for(k = 1;k < meta_atoms + 1;k++){ // *** temp_atomname1[0] = ""
      if( strstr( temp_atomname1[k], "TV") == NULL && strstr( temp_atomname1[k], "Tv") == NULL && temp_atomname1[k] != "" ){
        sprintf(atomname_TVexcluded[i], "%s", temp_atomname1[k]);
        i++;
      }/*
      printf("temp_atomname1 %d : %s\n",k,temp_atomname1[k]);
      printf("atomname_TVexcluded %d : %s\n",i,temp_atomname1[i]);*/
    }
    i = 0;
  totalTVnum = TVnum;
  totalEQinfonum++; // added 9424
  if( crystal_surface_mode == 1 && totalTVnum < 3 ) {
    printf("TV doesn't exist! Invalid EQ file...\n");
    return -1;
  }
  
  ////////////// element counting section <<begin>> (9821 moved)
  char totalelementname[atoms+1][N], psfstr[atoms+1][N];
  char temporary[atoms+1][N], atomname1[atoms+1][N], atomname2[atoms+1][atoms+1][N], tempelementname[N], tempelementname_for_count[N];
  for(m = 0;m < atoms+1 ;m++){
    forDOSinp_Species_order_list[m] = m;
  }
        
//  printf("\r  runing... (element counting)      ");
    sprintf(totalelementname[0],"%sX",atomname_TVexcluded[0]); ///// !!! element names of "totalelementname" array is "X " <! one space added !>
    totalelementnamenum++;
    for(j = 1;j < atoms;j++){ /// ((for all atoms))
      for(k = 0;k < totalelementnamenum ;k++){
        elementnamelength = strlen(totalelementname[k]);
        sprintf(tempelementname_for_count,"%.*s",elementnamelength-1,totalelementname[k]);
        if( strcmp(atomname_TVexcluded[j],tempelementname_for_count) == 0 ){ /// if some totalelementname[k] matches to temp_atomname1[j]
          addnewatom = 0;
        }
      }
      for(k = 0;k < totalelementnamenum ;k++){
        for(m = 0;m < totalelementnamenum ;m++){ /// ((for all element names))
          if( strcmp(totalelementname[m],totalelementname[k]) == 0 && m != k ){ /// if some totalelementname[k] overlup
            addnewatom = 0;
          }
        }
      }
      if( addnewatom == 1 ){
        sprintf(totalelementname[totalelementnamenum], "%sX", atomname_TVexcluded[j]);
        totalelementnamenum++;
      }
      addnewatom = 1; /// intitialize "addnewatom"
    }
    for(j = 0;j < totalelementnamenum;j++){
      for(k = j+1;k < totalelementnamenum ;k++){
        if( Zvalue(totalelementname[k]) < Zvalue(totalelementname[j]) ){
          sprintf(tempelementname, "%s", totalelementname[k]);
          sprintf(totalelementname[k], "%s", totalelementname[j]);
          sprintf(totalelementname[j], "%s", tempelementname);
        }
      }
    }
    tempnum = 0;
    for(k = 0;k < 200;k++){
      for(j = 0;j < atoms;j++){
        if( Zvalue(temp_atomname2[j]) == k ){
          if( Zvalue(temp_atomname2[j]) <= 10 ){
            sprintf(psfstr[tempnum], "  %d    %d    %s-pgpb", tempnum+1, Zvalue(temp_atomname2[j]), temp_atomname1[j] );
          }else if( Zvalue(temp_atomname2[j]) <= 18 ){
            sprintf(psfstr[tempnum], "  %d    %d    %s-pepb", tempnum+1, Zvalue(temp_atomname2[j]), temp_atomname1[j] );
          }else{
            sprintf(psfstr[tempnum], "  %d    %d    %s-pepbr", tempnum+1, Zvalue(temp_atomname2[j]), temp_atomname1[j] );
          }
          tempnum++;
          break;
        }
      }
    }
    tempnum = 0;
    for(k = 0;k < atoms+1;k++){
      for(j = 0;j < atoms+1;j++){
        if( Zvalue(temp_atomname2[k]) == Zvalue(temp_atomname2[j]) && forDOSinp_Species_order_list[k] > tempnum && forDOSinp_Species_order_list[j] > tempnum ){
          forDOSinp_Species_order_list[j] = tempnum;
          for(m = 0;m < atoms+1;m++){
            if( Zvalue(temp_atomname2[j]) == Zvalue(temp_atomname2[m]) && forDOSinp_Species_order_list[m] > tempnum ){
              forDOSinp_Species_order_list[m] = tempnum;
            }
          }
          tempnum++;
        }
      }
    }/*
    for(m = 1;m < atoms+1;m++){
        printf("%d\n",forDOSinp_Species_order_list[m]);
    }*/
  ////////////// element counting section <<end>> 
  
  printf("  JOBNAME : %s\n   atoms : %d, elements : %d, total mass : %.2lf", JOBNAME, atoms, totalelementnamenum, unit_cell_mass);
  if( logfile_mode == 0) {
    printf(",  the number of EQ : %d", totalEQnum);
  }
  printf("\n");
  fclose(fp);

  
  if( logfile_mode == 1 ){ //////######//////###### Path mode (BEGIN) //////######//////######
    char log_temporary[atoms+1][N];
  
    int IRCnum = 1,IRCing = 0;
    int OPTnum = 0,OPTing = 0;
    int FREQnum = 0,FREQing = 0;
    int Thermochemistry = 0,Eneprofile = 0;
    int AppTSnum = 0,LUPITRnum = -1,NODEnum = 0,latestNODEnum = 0;
    int FORWARD = 0,BACKWARD = 0;
    int MINoptimization = -1; // added 9430
    
  //  int G_TS[N], G_FOR[N], G_BACK[N];
    double G_TSval, G_FORval, G_BACKval;
    double temp_node_coordinate_x[L][atoms+totalTVnum], temp_node_coordinate_y[L][atoms+totalTVnum], temp_node_coordinate_z[L][atoms+totalTVnum], temp_node_energy[L]; // added 9430
    char temp_node_atomname[L][atoms+totalTVnum][N]; // added 9430
    char tempfilename[N];
    char max[N];
    char EXfilename[N];
    char forluptofile[K],echo[N],comstr[N];
    FILE *output2,*newcom,*oldcom,*IRCoutfile,*tempEQfile;
    char oldcomfilename[K];
    char newcomfilename[K];
    char infiles[N];
    char IRCoutfilename[K];

    fp = fopen( original_EQlist_fname,"r"); // open the original file to read (it is NOT EQlist !!)
    if(fp == NULL) { // failed to open, return NULL
      printf("%s.log doesn't exist!\n", JOBNAME);
      return -1;
    }
    tempEQfile = fopen( "temptemp_eq_EQ_list.log", "w"); // open the original file to read (it is NOT EQlist !!)

  // print date (BEGIN)
    char dateout[64];
    char datefilename[N];
    char tempdatenum[N],datenum[64];
    char optionsimbols[30];
    FILE *datefile;
    sprintf(datefilename, "date_%s.txt", JOBNAME);
    sprintf(optionsimbols, "+%%Y%%m%%d_%%H%%M");
    sprintf(dateout, "echo -n `date %s` > %s", optionsimbols, datefilename);
    system(dateout);
      datefile = fopen(datefilename,"r");
      fgets(tempdatenum, N, datefile);
      int datelength;
      datelength = strlen(datenum);
      sprintf(datenum, "%.*s", datelength - 1, tempdatenum);
      fclose(datefile);
    sprintf(dateout, "echo -n `rm date*.txt`");
    system(dateout);
  // print date (END)
  
    sprintf(filename, "EQCS_%s_EGLOG.xyz", JOBNAME);
    output = fopen(filename,"w"); //open the output file

    fprintf(output,"*************************************************************************\n");
    fprintf(output,"  Energy and Geometry profiles log\n");
    fprintf(output,"      \"%s\" (EQCS ver %s)\n", JOBNAME,software_version);
    fprintf(output,"*************************************************************************\n");
    fprintf(output,"                                            produced by H.N. (%s)\n\n",developed_date);

      while(fgets(str, N, fp) != NULL){
      
        if(strstr( str , "Geometry of AppTS" ) != NULL && IRCing == 0){ //"IRC..." found  (210708 modified)
          if(LUPITRnum > 0){ // first AppTS appeared (8827added)
            char LUPOUTname[K];
            FILE *LUPOUTout;
            sprintf(LUPOUTname, "%s_LUPOUTt.log", JOBNAME);
            LUPOUTout = fopen( LUPOUTname,"r");
            if(LUPOUTout != NULL) { //if LUPOUTt.log exists (8829refined) (9501revised)
              sprintf(forluptofile, "grep NODE -A%d %s_LUPOUTt.log > %s_forLUP_%s.log > 000_000_nu00", atoms, JOBNAME, JOBNAME, datenum); // output forLUP file (8827added)
              char *luppath = forluptofile;
              system(luppath);
              fclose(LUPOUTout);
            }
            
            // make new .com file for further LUP (8828added)
            sprintf(newcomfilename, "%s_nextLUP_%s.com", JOBNAME, datenum);
            sprintf(oldcomfilename, "%s.com", JOBNAME);
            oldcom = fopen(oldcomfilename,"r");
            if(oldcom != NULL) { // if LUPOUTt.log exists (8829refined) (9501revised)
              newcom = fopen(newcomfilename,"w");
              sprintf(infiles, "%%infile=%s_forLUP_%s", JOBNAME, datenum);
                fprintf(newcom, "%s\n", infiles);
                fgets(comstr, N, oldcom); // skip first sentence (the "%infile" line)
              while(fgets(comstr, N, oldcom) != NULL){
                fprintf(newcom, "%s", comstr);
              }
              fclose(oldcom);
              fclose(newcom);
            }
            //LUPITRnum = 0;
          }
          fprintf(output, " <<BEGIN>> Geometry of AppTS %d <<BEGIN>>\n=========================================================================\n=========================================================================\n",AppTSnum);
        }
        
        if(strstr( str , "INITIAL STRUCTURE" ) != NULL && IRCing == 0){ //"IRC..." found
          fprintf(output, "\n=========================================================================\n%d\nTS (%s)\n", meta_atoms, original_EQlist_fname); // (210708 modified)
          for(k = 0;k < meta_atoms;k++){ // (210708 added)
            fgets(str, N, fp);
            fprintf(output, "%s", str);
          }
          fprintf(output, "\n");
          IRCing = 1; //activate IRCing
        }
        if(strstr( str , "IRC FOLLOWING (FORWARD) STARTING FROM FIRST-ORDER SADDLE" ) != NULL){ //"FORWARD" begin
          FORWARD = 1; //activate FORWARD
        }
        if(strstr( str , "IRC FOLLOWING (BACKWARD) STARTING FROM FIRST-ORDER SADDLE" ) != NULL){ //"BACKWARD" begin
          BACKWARD = 1; //activate BACKWARD
          FORWARD = 0; //deactivate FORWARD
        }
        if(strstr( str , "LUP-path optimization" ) != NULL){ //count LUPITRnum (8827added)
          if( LUPITRnum < 0 ){
            LUPITRnum = 0;
          }
          latestNODEnum = NODEnum;
          NODEnum = 0;
          LUPITRnum++;
        }
        
        if(strstr( str , "# NODE" ) != NULL){ //count NODEnum (8c10added)
          for(k = 0;k < meta_atoms;k++){
            fgets(str, N, fp);
            sscanf(str, "%s    %lf    %lf    %lf\n", temp_node_atomname[NODEnum][k], &temp_node_coordinate_x[NODEnum][k], &temp_node_coordinate_y[NODEnum][k], &temp_node_coordinate_z[NODEnum][k]);
          }
          fgets(str, N, fp);
          sscanf(str, "ENERGY    = %lf\n", &temp_node_energy[NODEnum]);
          NODEnum++;
        }
        
        if(strstr( str , "Maximum number of iteration was exceeded" ) != NULL){ //if Max ITR exceeded
          fprintf(output, "\n%d\n", meta_atoms);
            if(FORWARD == 1){
              fprintf(output, "Max ITR was exceeded (# AppTS %d :FORWARD)\n", AppTSnum);
            }else{
              if(BACKWARD == 1){
                fprintf(output, "Max ITR was exceeded (# AppTS %d :BACKWARD)\n", AppTSnum);
              }else{
                fprintf(output, "Max ITR was exceeded (# AppTS %d)\n", AppTSnum);
              }
            }
            for(k = 0;k < meta_atoms;k++){
              fprintf(output, "H     0.000000000000     0.000000000000     0.000000000000\n");
            }
          fprintf(output, "\n <<Maximum number of iteration was exceeded (# AppTS %d)>>\n\n=========================================================================\n",AppTSnum);
          fprintf(output, "=========================================================================\n");
          fprintf(output, "\n\n\n");
          FILE *max; //to other file(s) (8827added)
          sprintf(EXfilename, "EQCS_%s_EGLOG_AppTS%d_exceeded.log", JOBNAME, AppTSnum);
          max = fopen(EXfilename, "w");
          fclose(max);
          IRCnum = IRCnum + 1; //IRCnum increase +1
          AppTSnum = AppTSnum + 1; //AppTSnum increase +1
          IRCing = 0; //deactivate IRCing
          FORWARD = 0; BACKWARD = 0; //initialize FORWARD and BACKWARD
        } //if Max ITR exceeded (end)
  
        if(strstr( str , "The number of peaks is too small..." ) != NULL){
          FILE *small; //to other file(s) (8827added)
          sprintf(EXfilename, "EQCS_%s_EGLOG_too_small_peaks.log", JOBNAME);
          small = fopen(EXfilename, "w");
          fclose(small);
        }
        if(strstr( str , "Couldn't find a lower energy point" ) != NULL){ // for IRC
          FILE *small; //to other file(s) (8827added)
          sprintf(EXfilename, "EQCS_%s_EGLOG_could_not_find_lowerEpoint.log", JOBNAME);
          small = fopen(EXfilename, "w");
          fclose(small);
        }
        
              //OPT part begin
              if(strstr( str , "MIN-optimization" ) != NULL){ //"OPT..." found
                MINoptimization++; // add MINoptimization
              }
              if(strstr( str , "OPTOPTOPT" ) != NULL){ //"OPT..." found
                OPTing = 1; //OPTing activate
              }
              if(strstr( str , "Optimized structure" ) != NULL && OPTing == 1){ //"Optimized structure..." found
                fprintf(output, "-------------------------------------------------------------------------\n");
                fprintf(output, "%d\n", meta_atoms);
                for(k = 0;k < meta_atoms + 1 ;k++){
                  fgets(log_temporary[k], N, fp);
                }
                if(FORWARD == 1){
                  fprintf(output, "FW (# AppTS %d) %s", AppTSnum,log_temporary[meta_atoms]);
                }else{
                  if(BACKWARD == 1){
                    fprintf(output, "BW (# AppTS %d) %s", AppTSnum,log_temporary[meta_atoms]);
                  }else if( MINoptimization < 0 ){
                    fprintf(output, "TS (# AppTS %d) %s", AppTSnum,log_temporary[meta_atoms]);
                  }else{
                    fprintf(output, "MINopt (AppEQ-%d)\n", MINoptimization);
                  }
                }
                for(k = 0;k < meta_atoms ;k++){
                  fprintf(output, "%s", log_temporary[k]);
                }
                OPTing = 0; //OPTing deactivate
                OPTnum = OPTnum + 1;
              }//OPT part (end)
                
              // if FREQ calculation was performed
              if(strstr( str , "Thermochemistry at" ) != NULL || strstr( str , "Thermochemistry (after" ) != NULL){ //"Thermochemistry" found (1629 modified)
                fprintf(output, "\n");
                  if(FORWARD == 1){
                    fprintf(output, "FW (# AppTS %d)", AppTSnum);
                  }else{
                    if(BACKWARD == 1){
                      fprintf(output, "BW (# AppTS %d)", AppTSnum);
                    }else if( MINoptimization < 0 ){
                      fprintf(output, "TS (# AppTS %d) %s", AppTSnum,log_temporary[meta_atoms]);
                    }else{
                      fprintf(output, "MINopt (AppEQ-%d)\n", MINoptimization);
                    }
                  }
                fprintf(output, "\n");
                for(k = 0;k < 15;k++){
                  fprintf(output, "%s", str);
                  if(strstr( str , "Free Energy" ) != NULL && FORWARD == 1){
                    sscanf(str,"Free Energy   = %lf\n", &G_FORval);
                  }else{
                    if(strstr( str , "Free Energy" ) != NULL && BACKWARD == 1){
                      sscanf(str,"Free Energy   = %lf\n", &G_BACKval);
                    }else if(strstr( str , "Free Energy" ) != NULL){
                      sscanf(str,"Free Energy   = %lf\n", &G_TSval);
                    }
                  }
                  fgets(str, N, fp);
                }
                fprintf(output, "\n");
                FREQnum = FREQnum + 1;
                FORWARD = 0; BACKWARD = 0; //initialize FORWARD and BACKWARD
              }// if FREQ calculation was performed (end)
  
              //Energy profile part begin
              if(strstr( str , "Energy profile along IRC" ) != NULL){ //Energy profile along IRC
                fgets(str, N, fp);                
//                fprintf(output, "\nEnergy profile along # IRC %d\n", IRCnum);
//                fprintf(output, "-------------------------------------------------------------------------\n");
//                fprintf(output, "\n");
                sprintf(IRCoutfilename, "EQCS_%s_EGLOG_AppTS%d_IRC.log", JOBNAME, AppTSnum);
                IRCoutfile = fopen(IRCoutfilename, "w");
                while(strstr( str , "Reverse" ) == NULL){ //until "Reverse" found
//                fprintf(output, "%s", str);
                  fprintf(IRCoutfile, "%s", str);
                  fgets(str, N, fp);
                }
                fclose(IRCoutfile);
                fprintf(output, "\n-------------------------------------------------------------------------\n");
                fprintf(output, "Energy profile along # IRC %d (# AppTS %d)\n", IRCnum, AppTSnum);
                fprintf(output, "IRC following along both forward and backward directions were finished\n");
                fprintf(output, "\n");
                fprintf(output, "FW =\t%17.12lf\nTS =\t%17.12lf\nBW =\t%17.12lf\n", G_FORval, G_TSval, G_BACKval);
                fprintf(output, "FW -> TS =\t%17.12lf\nBW -> TS =\t%17.12lf\n", 2625.49962*(G_TSval-G_FORval),2625.49962*(G_TSval-G_BACKval));
                fprintf(output, "\n");
                fprintf(output, " <<END>> IRC profiles of # AppTS %d <<END>>\n=========================================================================\n",AppTSnum);
                fprintf(output, "=========================================================================\n");
                fprintf(output, "\n\n\n\n");
                
                sprintf(tempfname, "EQCS_%s_EGLOG_FreeEnergy.log", JOBNAME);
                output2 = fopen(tempfname,"a+"); //open the output file (8c15 fixed)
                fprintf(output2, "Energy profile along # IRC %d (# AppTS %d)\n", IRCnum, AppTSnum);
                fprintf(output2, "-------------------------------------------------------------------------\n");
                fprintf(output2, "FW =\t%17.12lf\nTS =\t%17.12lf\nBW =\t%17.12lf\n", G_FORval, G_TSval, G_BACKval);
                fprintf(output2, "FW -> TS =\t%17.12lf\nBW -> TS =\t%17.12lf\n", 2625.49962*(G_TSval-G_FORval),2625.49962*(G_TSval-G_BACKval));
                fprintf(output2, "=========================================================================\n\n");
                fclose(output2);
                
                IRCnum = IRCnum + 1; //IRCnum increase +1
                AppTSnum = AppTSnum + 1; //AppTSnum increase +1
                IRCing = 0; //deactivate IRCing
                FORWARD = 0; BACKWARD = 0; //initialize FORWARD and BACKWARD
              }//Energy profile part (end)
  
        if(strstr( str , "termination" ) != NULL){
          fprintf(output, "\n");
          fprintf(output,"*************************************************************************\n");
          for(k = 0;k < 5;k++){
            fprintf(output, "%s", str);
            fgets(str, N, fp);
          }
          fprintf(output, "%s", str);
          fprintf(output,"*************************************************************************");
          break;
        }
  
      } //end of "while"
    
    FILE *LUPITRout,*NODEnumout;
    sprintf(tempfilename, "EQCS_%s_EGLOG_LUPITR_%d.log", JOBNAME, LUPITRnum);
    LUPITRout = fopen(tempfilename,"w");
    fclose(LUPITRout);
    
    sprintf(tempfilename, "EQCS_%s_EGLOG_NODE_%d.log", JOBNAME, latestNODEnum);
    NODEnumout = fopen(tempfilename,"w");
    fclose(NODEnumout);
    
    // make tempEQlist (begin)
    fprintf(tempEQfile, "List of Equilibrium Structures\n\n# Geometry of EQ 0, SYMMETRY = C1  \n");
    for(k = 0;k < meta_atoms;k++){
      fprintf(tempEQfile, "%-2s    %17.12lf    %17.12lf    %17.12lf\n", temp_node_atomname[0][k], temp_node_coordinate_x[0][k], temp_node_coordinate_y[0][k], temp_node_coordinate_z[0][k]);
    }
    fprintf(tempEQfile, "Energy    = %17.12lf\nSpin(**2) =    0.000000000000\nZPVE      =    0.000000000000\nNormal mode eigenvalues : nmode = %d\n", temp_node_energy[0], (meta_atoms * 3) - 6 );
    tempnum = ((meta_atoms * 3) - 6) / 5 ;
    for(k = 0;k < tempnum;k++){
      fprintf(tempEQfile, "  0.300000000   0.300000000   0.300000000   0.300000000   0.300000000\n");
    }
    tempnum = ((meta_atoms * 3) - 6) % 5 ;
    if( tempnum > 0){
      fprintf(tempEQfile, "  0.300000000");
      for(k = 1;k < tempnum;k++){
        fprintf(tempEQfile, "   0.300000000");
      }
      fprintf(tempEQfile, "\n");
    }
    
    fprintf(tempEQfile, "\n# Geometry of EQ 1, SYMMETRY = C1  \n");
    for(k = 0;k < meta_atoms;k++){
      fprintf(tempEQfile, "%-2s    %17.12lf    %17.12lf    %17.12lf\n", temp_node_atomname[latestNODEnum][k], temp_node_coordinate_x[latestNODEnum][k], temp_node_coordinate_y[latestNODEnum][k], temp_node_coordinate_z[latestNODEnum][k]);
    }
    fprintf(tempEQfile, "Energy    = %17.12lf\nSpin(**2) =    0.000000000000\nZPVE      =    0.000000000000\nNormal mode eigenvalues : nmode = %d\n", temp_node_energy[latestNODEnum], (meta_atoms * 3) - 6 );
    tempnum = ((meta_atoms * 3) - 6) / 5 ;
    for(k = 0;k < tempnum;k++){
      fprintf(tempEQfile, "  0.300000000   0.300000000   0.300000000   0.300000000   0.300000000\n");
    }
    tempnum = ((meta_atoms * 3) - 6) % 5 ;
    if( tempnum > 0){
      fprintf(tempEQfile, "  0.300000000");
      for(k = 1;k < tempnum;k++){
        fprintf(tempEQfile, "   0.300000000");
      }
      fprintf(tempEQfile, "\n");
    }
    
    // check NODECONNECTION
    FILE *NODECONNECTION;
    if( no_babel_mode == 1 ){
      sprintf(command, "%s temptemp_eq -c -n > 000_000_nu00", argv[0]);
      system(command);
      sprintf(tempfilename, "EQCS_temptemp_eq_cLOG.log");
      NODECONNECTION = fopen(tempfilename,"r");
      while(fgets(str, N, NODECONNECTION) != NULL){
        if(strstr( str , "the number of composition" ) != NULL){
          sscanf(str, "   the number of composition     : %d\n", &tempnum);
          break;
        }
      }
      fclose(NODECONNECTION);
      if(tempnum == 2){
        sprintf(tempfilename, "EQCS_%s_EGLOG_diff_type.log", JOBNAME);
      }else{
        sprintf(tempfilename, "EQCS_%s_EGLOG_same_type.log", JOBNAME);
      }
      NODECONNECTION = fopen(tempfilename,"w");
      fclose(NODECONNECTION);
    }else if( no_babel_mode == 0 ){
      // sprintf(command, "%s temptemp_eq -c > 000_000_nu00", argv[0]);
      sprintf(command, "%s temptemp_eq -c > 000_000_nu00", argv[0]);
      system(command);
      sprintf(tempfilename, "EQCS_temptemp_eq_cLOG.log");
      NODECONNECTION = fopen(tempfilename,"r");
      while(fgets(str, N, NODECONNECTION) != NULL){
        if(strstr( str , "the number of InChI sequence" ) != NULL){
          sscanf(str, "   the number of InChI sequence  : %d\n", &tempnum2);
          fgets(str, N, NODECONNECTION);
          sscanf(str, "   the number of composition     : %d\n", &tempnum);
          break;
        }
      }
      fclose(NODECONNECTION);
      if(tempnum == 2){
        if(tempnum2 == 2){
          sprintf(tempfilename, "EQCS_%s_EGLOG_diff_type_diff_seq.log", JOBNAME);
        }else{
          sprintf(tempfilename, "EQCS_%s_EGLOG_diff_type_same_seq.log", JOBNAME);
        }
      }else{
        if(tempnum2 == 2){
          sprintf(tempfilename, "EQCS_%s_EGLOG_same_type_diff_seq.log", JOBNAME);
        }else{
          sprintf(tempfilename, "EQCS_%s_EGLOG_same_type_same_seq.log", JOBNAME);
        }
      }
      NODECONNECTION = fopen(tempfilename,"w");
      fclose(NODECONNECTION);

      NODECONNECTION = fopen(tempfilename,"w");
      fclose(NODECONNECTION);
    }
    // make tempEQlist (end)
    fclose(tempEQfile); //close the tempEQlist file
    fclose(output); //close the output file
    fclose(fp); //close the original file
    // delete intermediate files
    remove("temptemp_eq_EQ_list.log");
    remove("EQCS_temptemp_eq_EQLOG.xyz"); // 1630 added
    remove("EQCS_temptemp_eq_EQlistLOG.xyz");
    remove("EQCS_temptemp_eq_cLOG.log");
    remove("EQCS_temptemp_eq_bLOG.log");
    //
    printf("\r  EG sorting was completed. ");
    printf("(EQCS ver %s)\n", software_version);
//    return 0;
    if( composition_sorting_option == 1 ){
      if(EQlist2_mode == 0 ) {
        sprintf( original_EQlist_fname,"%s_EQ_list.log",JOBNAME);
        fp = fopen( original_EQlist_fname,"r"); // open the original file to read
        if(fp == NULL) { // failed to open, return NULL
          sprintf( original_EQlist_fname,"%s_EQ_list_2.log",JOBNAME);
          fp = fopen( original_EQlist_fname,"r"); // open the original file to read
          if(fp == NULL) { // failed to open, return NULL
            printf("%s_EQ_list.log and list_2.log doesn't exist!\n", JOBNAME);
            return -1;
          }else{
            EQlist2_mode = 1;
          }
        }
      }else if(EQlist2_mode == 1) {
        sprintf( original_EQlist_fname,"%s_EQ_list_2.log",JOBNAME);
        fp = fopen( original_EQlist_fname,"r"); // open the original file to read
        if(fp == NULL) { // failed to open, return NULL
          printf("%s_EQ_list_2.log doesn't exist!\n", JOBNAME);
          return -1;
        }
      }
      while( fgets(str, N, fp) != NULL ){
        if( strstr(str, "# Geometry of EQ" ) != NULL ){
          totalEQnum++;
        }
      }
      fclose(fp);
    }else{
      return 0;
    }
  } //////######//////###### Path mode (END) //////######//////######
  
  tempnum = 1;

  fp = fopen( original_EQlist_fname,"r"); //open the original file to read
  if(fp == NULL) { //failed to open, return NULL
    printf("failed to open %s!\n", original_EQlist_fname);
    return -1;
  }
  
  sprintf(xyzfilename, "EQCS_%s_EQLOG.xyz", JOBNAME);
  if(EQlist2_mode == 1) {
    sprintf(xyzfilename, "EQCS_%s_EQLOG2.xyz", JOBNAME);
  }
  output = fopen(xyzfilename,"w"); //open the output file
  if(output == NULL) { //failed to open, return NULL
    printf( "failed to open %s!\n", xyzfilename);
    return -1;
  }
  
  if(energy_sort_mode == 1 && EQlist2_mode == 0){
    sprintf(elogfilename, "EQCS_%s_eLOG.log", JOBNAME);
    energy_sort = fopen(elogfilename,"w"); //open the output file
    if(energy_sort == NULL) { //failed to open, return NULL
      printf( "failed to open %s!\n", elogfilename);
      return -1;
    }
    sprintf(exyzfilename, "EQCS_%s_eLOG.xyz", JOBNAME);
    energy_sort_xyz = fopen(exyzfilename,"w"); //open the output file
    if(energy_sort_xyz == NULL) { //failed to open, return NULL
      printf( "failed to open %s!\n", exyzfilename);
      return -1;
    }
  }
  
  //////////// 9615 added (for Cytoscape)
  if( totalEQnum > 0 && cytoscapelistgen_option == 1 ){
    sprintf(cytofilename, "EQCS_%s_cytoEQLOG.csv", JOBNAME);
    cytoscapeEQlist = fopen(cytofilename,"w"); //open the output file
    if(cytoscapeEQlist == NULL) { //failed to open, return NULL
      printf( "failed to open %s!\n", cytofilename);
      return -1;
    }
  }
  
  if(density_sort_mode == 1 && crystal_surface_mode == 1){
    sprintf(dlogfilename, "EQCS_%s_dLOG.log", JOBNAME);
    density_sort = fopen(dlogfilename,"w"); //open the output file
    if(density_sort == NULL) { //failed to open, return NULL
      printf( "failed to open %s!\n", dlogfilename);
      return -1;
    }
  }
  
  if( crystal_surface_mode == 1 ){
    sprintf(scaled_EQlist_filename, "EQCS_%s_scaled_EQ_list.log", JOBNAME);
    scaledEQlist = fopen(scaled_EQlist_filename,"w"); //open the output file
    fprintf(scaledEQlist, "List of Equilibrium Structures\n\n");
    if(scaledEQlist == NULL) { //failed to open, return NULL
      printf( "failed to open %s!\n", scaled_EQlist_filename);
      return -1;
    }
    sprintf(cell_VDlist_filename, "EQCS_%s_Volume_Density_list.log", JOBNAME);
    cellVDlist = fopen(cell_VDlist_filename,"w"); //open the output file
    if(cellVDlist == NULL) { //failed to open, return NULL
      printf( "failed to open %s!\n", cell_VDlist_filename);
      return -1;
    }
    fprintf(cellVDlist, "List of Volume and Density (%s)\nNo.",JOBNAME);
    for(k = 0;k < totalEQ_digitnum;k++){
      fprintf(cellVDlist, " ");
    }
    fprintf(cellVDlist, " : volume per cell (A^3) : density (g/cm^3)\n");
    
  }
  
  // for super cell generating (9613 added)
  if( supercellgen_option == 1 && crystal_surface_mode == 1 ){
    sprintf(spfilename, "EQCS_%s_SPcell_EQ_list.log", JOBNAME);
    supercellEQlist = fopen(spfilename,"w"); //open the output file
    if(supercellEQlist == NULL) { //failed to open, return NULL
      printf( "failed to open %s!\n", spfilename);
      return -1;
    }
    sprintf(spfilename, "EQCS_%s_SPcell.xyz", JOBNAME);
    supercellxyz = fopen(spfilename,"w"); //open the output file
    if(supercellxyz == NULL) { //failed to open, return NULL
      printf( "failed to open %s!\n", spfilename);
      return -1;
    }
  }
/*  作りかけ(9728)
  if( rearranged_EQlist_output_option == 1 ){
    sprintf(rearranged_EQlist_filename, "EQCS_%s_arr_EQ_list.log", JOBNAME);
    rearrangedEQlist = fopen(rearranged_EQlist_filename,"w"); //open the output file
    if(rearrangedEQlist == NULL) { //failed to open, return NULL
      printf( "failed to open %s!\n", rearranged_EQlist_filename);
      return -1;
    }
  }*/
  
//  printf("%d\n",totalEQinfonum);
/*
  fprintf(output,"*************************************************************************\n");
  fprintf(output,"  EQ list sorting log of\n");
  fprintf(output,"      \"%s\" (ver %s)\n", JOBNAME,software_version);
  fprintf(output,"*************************************************************************\n");
  fprintf(output,"                                        produced by H.Nabata (%s)\n\n",developed_date);
*/
/// variables for EQ list up (arrays required as many as totalEQnum)
  double Forced_Energy[totalEQnum], Bare_Energy[totalEQnum], sorted_Bare_Energy[totalEQnum], Zero_Energy[totalEQnum];
  double Electron_Energy[totalEQnum], zero_Gibbs_Energy[totalEQnum], RT_Gibbs_Energy[totalEQnum], Boltzmann_Distribution[totalEQnum]; // for EQ_list_2 (9411 added)
  double cell_length_a[totalEQnum], cell_length_b[totalEQnum], cell_length_c[totalEQnum], alpha_yz[totalEQnum], beta_xz[totalEQnum], gamma_xy[totalEQnum], v_for_cell[totalEQnum]; // (9715 added)
  char atomname[totalEQnum][atoms+totalTVnum][N], temporary_TVincluded[atoms+totalTVnum+1][N], arr_temporary_TVincluded[atoms+totalTVnum+1][N], EQinfo[totalEQnum][totalEQinfonum][N], EQinfo_temp_char1[N], EQinfo_temp_char2[N], topline[N];
  char elemetname_list[atoms][N]; // added 9821
  int shortest_TVnum[totalEQnum], spEQnum = 0; // added 9613
  double temp_coordinate_x[totalEQnum][atoms+totalTVnum], temp_coordinate_y[totalEQnum][atoms+totalTVnum], temp_coordinate_z[totalEQnum][atoms+totalTVnum], cell_Vol[totalEQnum], cell_density[totalEQnum], sorted_cell_density[totalEQnum], temp_shortestTVlength = 100000.0;
//  double arr_coordinate_x[totalEQnum][atoms+totalTVnum];
  int oldEQnum[totalEQnum], Energy_rank_EQnum[totalEQnum], Density_rank_EQnum[totalEQnum];
    for(k = 0;k < totalEQnum ;k++){ /// intialize "bonding_groups"
      Energy_rank_EQnum[k] = k;
      Density_rank_EQnum[k] = k;
    }
    
  while(fgets(str, N, fp) != NULL){
    if( strstr( str , "#" ) != NULL ){
      if( gjfgen_option == 1 ){ // 9613 added
        sprintf(gjffilename, "EQCS_%s_EQ%d.gjf", JOBNAME, spEQnum);
        gjf = fopen(gjffilename,"w"); //open the output file
        if(gjf == NULL) { //failed to open, return NULL
          printf( "failed to open %s!\n", gjffilename);
          return -1;
        }
        fprintf(gjf, "%%chk=%s_EQ%d.chk\n# pbepbe/6-31g/auto\n\n%s_EQ%d (%d %d %d)\n\n0 1\n", JOBNAME, tempEQnum, JOBNAME, tempEQnum, xyzframescale_x, xyzframescale_y, xyzframescale_z);
      }
      if( cifgen_option == 1 && crystal_surface_mode == 1 ){ // 9715 added (FOR CIF FILES)
        sprintf(ciffilename, "EQCS_%s_EQ%d.cif", JOBNAME, spEQnum);
        cif = fopen(ciffilename,"w"); //open the output file
        if(cif == NULL) { //failed to open, return NULL
          printf( "failed to open %s!\n", ciffilename);
          return -1;
        }
        if( spEQnum == 0 ){
          fprintf(cif, "data_%s_EQ%d\n", JOBNAME, spEQnum);
          sprintf(ciflistfilename, "EQCS_%s_EQ_list.cif", JOBNAME ); // 9928 added (FOR CIF_LIST )
          ciflist = fopen(ciflistfilename,"w"); //open the output file
          if(ciflist == NULL) { //failed to open, return NULL
            printf( "failed to open %s!\n", ciflistfilename);
            return -1;
          }
          fprintf(ciflist, "data_%s_EQ%d\n", JOBNAME, spEQnum);
        }else{
          fprintf(ciflist, "data_%s_EQ%d\n", JOBNAME, spEQnum);
        }
      }
      if( forDOSinpgen_option == 1 && crystal_surface_mode == 1 ){ // ( DOS inp #1 BEGIN)
        sprintf(forDOSinpfilename, "DOS_%s_EQ%d.inp", JOBNAME, tempEQnum);
        forDOSinp = fopen(forDOSinpfilename,"w"); //open the output file
        if(forDOSinp == NULL) { //failed to open, return NULL
          printf( "failed to open %s!\n", forDOSinpfilename);
          return -1;
        }
        sprintf(forDOSshfilename, "DOS_%s_EQ%d.sh", JOBNAME, tempEQnum);
        forDOSsh = fopen(forDOSshfilename,"w"); //open the output file
        if(forDOSsh == NULL) { //failed to open, return NULL
          printf( "failed to open %s!\n", forDOSshfilename);
          return -1;
        }
        fprintf(forDOSinp, "SystemName        DOS_%s_EQ%d          # Descriptive name of the system\nSystemLabel      DOS_%s_EQ%d          # Short name for naming files\n\nNumberOfAtoms            %d           # Number of atoms\nNumberOfSpecies          %d           # Number of species\n\n%%block Chemical_Species_Label\n", JOBNAME, tempEQnum, JOBNAME, tempEQnum, atoms, totalelementnamenum );
        for(j = 0;j < totalelementnamenum;j++){
          fprintf(forDOSinp, "%s\n",psfstr[j]);
        }
        fprintf(forDOSinp, "%%endblock Chemical_Species_Label\n\n\n# Lattice, coordinates, k-sampling\n\nAtomicCoordinatesFormat     Ang\n%%block AtomicCoordinatesAndAtomicSpecies\n");
      } // ( DOS inp #1 END)
      
      sprintf(topline, "%s",str);
      lntrim(topline);
      if( EQlist2_mode == 1 ){
        sscanf(str, "# Geometry of EQ %s (Old number = %d), SYMMETRY = %s\n", EQinfo_temp_char1, &oldEQnum[tempEQnum], EQinfo_temp_char2);
      }
      if( crystal_surface_mode == 0 ){
        fprintf(output, "%d\n", atoms);
      }else if( crystal_surface_mode == 1 ){
        fprintf(output, "%d\n", atoms * (framescale_x+1) * (framescale_y+1) * (framescale_z+1) );
        fprintf(scaledEQlist, "%s",str);
        if( supercellgen_option == 1 ){
          fprintf(supercellxyz, "%d\n", atoms * 2 );
          fprintf(supercellEQlist, "%s", str );
          sprintf(spfilename, "EQCS_%s_EQ%d.gjf", JOBNAME, spEQnum);
        }
      }
      for(k = 0;k < meta_atoms;k++){
        fgets(temporary_TVincluded[k], N, fp);
      }
      if( rearranged_EQlist_output_option == 1 ){ // added 9728
        for(m = 0;m < meta_atoms;m++){
          sprintf( arr_temporary_TVincluded[m], "%s", temporary_TVincluded[m]);
        }
        tempnum = 0;
        for(m = 0;m < meta_atoms;m++){
          if( strstr( temporary_TVincluded[m], "TV") == NULL && strstr( temporary_TVincluded[m], "Tv") == NULL ){
            sprintf( arr_temporary_TVincluded[m], "%s", temporary_TVincluded[rearranged_atom_order[m+1]-1+tempnum]);
            //printf(  "%s", temporary_TVincluded[rearranged_atom_order[m+1]-1]);
          }else{
            tempnum++;
          }
        }
        tempnum = 0;
        for(m = 0;m < meta_atoms;m++){
          if( strstr( arr_temporary_TVincluded[m], "TV") == NULL && strstr( arr_temporary_TVincluded[m], "Tv") == NULL ){
            sprintf( temporary_TVincluded[m], "%s", arr_temporary_TVincluded[m]);
          }else{
            tempnum++;
          }
        }
      }
      for(k = 0;k < meta_atoms;k++){
        sscanf(temporary_TVincluded[k], "%s    %lf    %lf    %lf\n", atomname[tempEQnum][k], &temp_coordinate_x[tempEQnum][k], &temp_coordinate_y[tempEQnum][k], &temp_coordinate_z[tempEQnum][k]);
      }
      fgets(EQinfo[tempEQnum][0], N, fp);
      if( strstr( EQinfo[tempEQnum][0] , "(" ) == NULL ){
        LUP_mode = 1;
      }
      if( EQlist2_mode == 0 && LUP_mode == 0 ){
        sscanf(EQinfo[tempEQnum][0], "Energy    = %lf ( %lf : lf)\n", &Forced_Energy[tempEQnum], &Bare_Energy[tempEQnum], &Zero_Energy[tempEQnum]);
        if( Bare_Energy[tempEQnum] > -0.1 ){
          Bare_Energy[tempEQnum] = Forced_Energy[tempEQnum];
        }
        fprintf(output, "# EQ %d ; %17.12lf (%17.12lf) ; %s\n", tempEQnum, Forced_Energy[tempEQnum], Bare_Energy[tempEQnum], topline);
        if( supercellgen_option == 1 ){
          fprintf(supercellxyz, "# EQ %d ; %17.12lf (%17.12lf) ; %s\n", tempEQnum, Forced_Energy[tempEQnum], Bare_Energy[tempEQnum], topline);
        }
/*        if( supercellgen_option == 1 ){
          fprintf(rearrangedEQlist, "%s", topline);
        }*/
      }else if( EQlist2_mode == 0 && LUP_mode == 1 ){
        sscanf(EQinfo[tempEQnum][0], "Energy    = %lf\n", &Bare_Energy[tempEQnum]);
        fprintf(output, "# EQ %d ; %17.12lf\n", tempEQnum, Bare_Energy[tempEQnum]);
      }else if( EQlist2_mode == 1 ){
        sscanf(EQinfo[tempEQnum][0], "Relat.Energy               = %lf kJ/mol\n", &Electron_Energy[tempEQnum]);
        fgets(EQinfo[tempEQnum][1], N, fp);
        sscanf(EQinfo[tempEQnum][1], "Relat.Energy   (   0.00 K) = %lf kJ/mol\n", &zero_Gibbs_Energy[tempEQnum]);
        fgets(EQinfo[tempEQnum][2], N, fp);
        sscanf(EQinfo[tempEQnum][2], "Relat.Energy   ( 298.15 K) = %lf kJ/mol\n", &RT_Gibbs_Energy[tempEQnum]);
        fgets(EQinfo[tempEQnum][3], N, fp);
        sscanf(EQinfo[tempEQnum][3], "Boltzmann Dist.( 298.15 K) = %lf %%\n", &Boltzmann_Distribution[tempEQnum]);
        fprintf(output, "# EQ %d (Old EQ %d) Rel.E=%.1f; Rel.G(0K)=%.1f; Rel.G(RT)=%.1f; (%.3f%%)\n", tempEQnum, oldEQnum[tempEQnum], Electron_Energy[tempEQnum], zero_Gibbs_Energy[tempEQnum], RT_Gibbs_Energy[tempEQnum], Boltzmann_Distribution[tempEQnum]);
        fgets(EQinfo[tempEQnum][4], N, fp);
        fgets(EQinfo[tempEQnum][5], N, fp);
      }
      if( crystal_surface_mode == 0 ){ // when crystal or surface mode is off
        for(k = 0;k < meta_atoms;k++){
          if( strstr( atomname[tempEQnum][k], "TV") == NULL && strstr( atomname[tempEQnum][k], "Tv") == NULL ){
            fprintf(output, "%s", temporary_TVincluded[k]);
            
            if( gjfgen_option == 1 ){
              fprintf(gjf, " %s", temporary_TVincluded[k]);
            }
          }
        }
      }else if( crystal_surface_mode == 1 ){ /////////// when crystal or surface mode is on (begin)
        if( EQlist2_mode == 0 ){
          for(k = 0;k < totalEQinfonum;k++){
            fgets(EQinfo[tempEQnum][k], N, fp);
          }
          sscanf(EQinfo[tempEQnum][totalEQinfonum-1], "Normal mode eigenvalues : nmode = %d\n", &nmodenum);
        }else if( EQlist2_mode == 1 ){
          sscanf(EQinfo[tempEQnum][5], "Harmonic frequencies : nmode = %d\n", &nmodenum);
        }
        if( TVextension_option == 1 ){ // TVextension (added 9a24)
          temp_coordinate_x[tempEQnum][TVatom_linenum[0]] = temp_coordinate_x[tempEQnum][TVatom_linenum[0]] + TVextension;
        }
        if( cifgen_option == 1 ){
          cell_length_a[tempEQnum] = pow(temp_coordinate_x[tempEQnum][TVatom_linenum[0]]*temp_coordinate_x[tempEQnum][TVatom_linenum[0]] + temp_coordinate_y[tempEQnum][TVatom_linenum[0]]*temp_coordinate_y[tempEQnum][TVatom_linenum[0]] + temp_coordinate_z[tempEQnum][TVatom_linenum[0]]*temp_coordinate_z[tempEQnum][TVatom_linenum[0]], 0.5);
          cell_length_b[tempEQnum] = pow(temp_coordinate_x[tempEQnum][TVatom_linenum[1]]*temp_coordinate_x[tempEQnum][TVatom_linenum[1]] + temp_coordinate_y[tempEQnum][TVatom_linenum[1]]*temp_coordinate_y[tempEQnum][TVatom_linenum[1]] + temp_coordinate_z[tempEQnum][TVatom_linenum[1]]*temp_coordinate_z[tempEQnum][TVatom_linenum[1]], 0.5);
          cell_length_c[tempEQnum] = pow(temp_coordinate_x[tempEQnum][TVatom_linenum[2]]*temp_coordinate_x[tempEQnum][TVatom_linenum[2]] + temp_coordinate_y[tempEQnum][TVatom_linenum[2]]*temp_coordinate_y[tempEQnum][TVatom_linenum[2]] + temp_coordinate_z[tempEQnum][TVatom_linenum[2]]*temp_coordinate_z[tempEQnum][TVatom_linenum[2]], 0.5);
          if( cifbohr_option == 1 ){
            fprintf(cif, "_cell_length_a %lf\n", cell_length_a[tempEQnum] * (xyzframescale_x + 1) / 1.8897259886); // in bohr (23829)
            fprintf(cif, "_cell_length_b %lf\n", cell_length_b[tempEQnum] * (xyzframescale_y + 1) / 1.8897259886);
            fprintf(cif, "_cell_length_c %lf\n", cell_length_c[tempEQnum] * (xyzframescale_y + 1) / 1.8897259886);
            fprintf(ciflist, "_cell_length_a %lf\n", cell_length_a[tempEQnum] * (xyzframescale_x + 1) / 1.8897259886);
            fprintf(ciflist, "_cell_length_b %lf\n", cell_length_b[tempEQnum] * (xyzframescale_y + 1) / 1.8897259886);
            fprintf(ciflist, "_cell_length_c %lf\n", cell_length_c[tempEQnum] * (xyzframescale_y + 1) / 1.8897259886);
          }else{
            fprintf(cif, "_cell_length_a %lf\n", cell_length_a[tempEQnum] * (xyzframescale_x + 1)); // in ang (23829)
            fprintf(cif, "_cell_length_b %lf\n", cell_length_b[tempEQnum] * (xyzframescale_y + 1));
            fprintf(cif, "_cell_length_c %lf\n", cell_length_c[tempEQnum] * (xyzframescale_y + 1));
            fprintf(ciflist, "_cell_length_a %lf\n", cell_length_a[tempEQnum] * (xyzframescale_x + 1));
            fprintf(ciflist, "_cell_length_b %lf\n", cell_length_b[tempEQnum] * (xyzframescale_y + 1));
            fprintf(ciflist, "_cell_length_c %lf\n", cell_length_c[tempEQnum] * (xyzframescale_y + 1));
          }
          alpha_yz[tempEQnum] = acos((temp_coordinate_x[tempEQnum][TVatom_linenum[2]]*temp_coordinate_x[tempEQnum][TVatom_linenum[1]] + temp_coordinate_y[tempEQnum][TVatom_linenum[2]]*temp_coordinate_y[tempEQnum][TVatom_linenum[1]] + temp_coordinate_z[tempEQnum][TVatom_linenum[2]]*temp_coordinate_z[tempEQnum][TVatom_linenum[1]])/(cell_length_b[tempEQnum] * cell_length_c[tempEQnum]));
          beta_xz[tempEQnum] = acos((temp_coordinate_x[tempEQnum][TVatom_linenum[0]]*temp_coordinate_x[tempEQnum][TVatom_linenum[2]] + temp_coordinate_y[tempEQnum][TVatom_linenum[0]]*temp_coordinate_y[tempEQnum][TVatom_linenum[2]] + temp_coordinate_z[tempEQnum][TVatom_linenum[0]]*temp_coordinate_z[tempEQnum][TVatom_linenum[2]])/(cell_length_a[tempEQnum] * cell_length_c[tempEQnum]));
          gamma_xy[tempEQnum] = acos((temp_coordinate_x[tempEQnum][TVatom_linenum[0]]*temp_coordinate_x[tempEQnum][TVatom_linenum[1]] + temp_coordinate_y[tempEQnum][TVatom_linenum[0]]*temp_coordinate_y[tempEQnum][TVatom_linenum[1]] + temp_coordinate_z[tempEQnum][TVatom_linenum[0]]*temp_coordinate_z[tempEQnum][TVatom_linenum[1]])/(cell_length_a[tempEQnum] * cell_length_b[tempEQnum]));
          fprintf(cif, "_cell_angle_alpha %lf\n_cell_angle_beta %lf\n_cell_angle_gamma %lf\n_chemical_formula_sum	%s_EQ%d\n\n_symmetry_space_group_name_H-M    \"P 1\"\n_symmetry_int_tables_number       1\n\nloop_\n	_symmetry_equiv_pos_as_xyz\n	'x, y, z'\n\nloop_\n	_atom_site_label\n	_atom_site_occupancy\n	_atom_site_fract_x\n	_atom_site_fract_y\n	_atom_site_fract_z\n	_atom_site_thermal_displace_type\n	_atom_site_B_iso_or_equiv\n	_atom_site_type_symbol\n", alpha_yz[tempEQnum] * 180.0 / (atan(1.0) * 4.0), beta_xz[tempEQnum] * 180.0 / (atan(1.0) * 4.0), gamma_xy[tempEQnum] * 180.0 / (atan(1.0) * 4.0), JOBNAME, spEQnum);
          fprintf(ciflist, "_cell_angle_alpha %lf\n_cell_angle_beta %lf\n_cell_angle_gamma %lf\n_chemical_formula_sum	%s_EQ%d\n\n_symmetry_space_group_name_H-M    \"P 1\"\n_symmetry_int_tables_number       1\n\nloop_\n	_symmetry_equiv_pos_as_xyz\n	'x, y, z'\n\nloop_\n	_atom_site_label\n	_atom_site_occupancy\n	_atom_site_fract_x\n	_atom_site_fract_y\n	_atom_site_fract_z\n	_atom_site_thermal_displace_type\n	_atom_site_B_iso_or_equiv\n	_atom_site_type_symbol\n", alpha_yz[tempEQnum] * 180.0 / (atan(1.0) * 4.0), beta_xz[tempEQnum] * 180.0 / (atan(1.0) * 4.0), gamma_xy[tempEQnum] * 180.0 / (atan(1.0) * 4.0), JOBNAME, spEQnum);
          v_for_cell[tempEQnum] = pow(1 - cos(alpha_yz[tempEQnum]) * cos(alpha_yz[tempEQnum]) - cos(beta_xz[tempEQnum]) * cos(beta_xz[tempEQnum]) - cos(gamma_xy[tempEQnum]) * cos(gamma_xy[tempEQnum]) + 2 * cos(alpha_yz[tempEQnum]) * cos(beta_xz[tempEQnum]) * cos(gamma_xy[tempEQnum]), 0.5);
        }

        l = 0;
        tempnum = 1; tempnum2 = 1;
        for(k = 0;k < meta_atoms;k++){
          if( strstr( atomname[tempEQnum][k], "TV") == NULL && strstr( atomname[tempEQnum][k], "Tv") == NULL ){
            fprintf(output, "%s",temporary_TVincluded[k]);
            fprintf(scaledEQlist, "%s",temporary_TVincluded[k]);
//            fprintf(spEQfile, "%s",temporary_TVincluded[k]);
            if( gjfgen_option == 1 ){
              fprintf(gjf, "%s", temporary_TVincluded[k]);
            }
            if( cifgen_option == 1 ){
              sprintf(temp_buff, "      %s%d  1.0000  %lf    %lf    %lf Biso 1.000 %s\n", atomname[tempEQnum][k], tempnum, temp_coordinate_x[tempEQnum][k] / cell_length_a[tempEQnum] - temp_coordinate_y[tempEQnum][k] * cos(gamma_xy[tempEQnum]) / (cell_length_a[tempEQnum] * sin(gamma_xy[tempEQnum])) + temp_coordinate_z[tempEQnum][k] * (cos(alpha_yz[tempEQnum]) * cos(gamma_xy[tempEQnum]) - cos(beta_xz[tempEQnum])) / (cell_length_a[tempEQnum] * v_for_cell[tempEQnum] * sin(gamma_xy[tempEQnum])) , temp_coordinate_y[tempEQnum][k] / (cell_length_b[tempEQnum] * sin(gamma_xy[tempEQnum])) + temp_coordinate_z[tempEQnum][k] * (cos(beta_xz[tempEQnum]) * cos(gamma_xy[tempEQnum]) - cos(alpha_yz[tempEQnum])) / (cell_length_b[tempEQnum] * v_for_cell[tempEQnum] * sin(gamma_xy[tempEQnum])) , temp_coordinate_z[tempEQnum][k] * sin(gamma_xy[tempEQnum]) / (cell_length_c[tempEQnum] * v_for_cell[tempEQnum]),  atomname[tempEQnum][k]);
              fprintf(cif, "%s", temp_buff);
              fprintf(ciflist, "%s", temp_buff);
              tempnum++;
            }
            if( forDOSinpgen_option == 1 ){ // ( DOS inp #2 BEGIN )
              fprintf(forDOSinp, "%17.12lf %17.12lf %17.12lf %d\n",temp_coordinate_x[tempEQnum][k], temp_coordinate_y[tempEQnum][k], temp_coordinate_z[tempEQnum][k], forDOSinp_Species_order_list[tempnum2] + 1);
              tempnum2++;
            } // ( DOS inp #2 END )
          }else{
            if( temp_shortestTVlength > temp_coordinate_x[tempEQnum][k]*temp_coordinate_x[tempEQnum][k] + temp_coordinate_y[tempEQnum][k]*temp_coordinate_y[tempEQnum][k] + temp_coordinate_z[tempEQnum][k]*temp_coordinate_z[tempEQnum][k] ){
              shortest_TVnum[tempEQnum] = l;
              temp_shortestTVlength = temp_coordinate_x[tempEQnum][k]*temp_coordinate_x[tempEQnum][k] + temp_coordinate_y[tempEQnum][k]*temp_coordinate_y[tempEQnum][k] + temp_coordinate_z[tempEQnum][k]*temp_coordinate_z[tempEQnum][k];
            }
            l++;
          }
        }
        if( cifgen_option == 1 ){
          fprintf(cif, "\n\n");
          fprintf(ciflist, "\n\n");
          fclose(cif);
        }
        temp_shortestTVlength = 100000.0; // initialize "temp_shortestTVlength"
        tempnum = 0; tempnum2 = 0;
          /* ( 0 ~ +framescale ) */
        for(k = 0;k <= framescale_x;k++){
          for(l = 0;l <= framescale_y;l++){
            for(m = 0;m <= framescale_z;m++){
              for(j = 0;j < meta_atoms;j++){
                if( j != TVatom_linenum[0] && j != TVatom_linenum[1] && j != TVatom_linenum[2] && k+l+m != 0 ){
                  double tempmeta_X = temp_coordinate_x[tempEQnum][j] + k*temp_coordinate_x[tempEQnum][TVatom_linenum[0]] + l*temp_coordinate_x[tempEQnum][TVatom_linenum[1]] + m*temp_coordinate_x[tempEQnum][TVatom_linenum[2]];
                  double tempmeta_Y = temp_coordinate_y[tempEQnum][j] + k*temp_coordinate_y[tempEQnum][TVatom_linenum[0]] + l*temp_coordinate_y[tempEQnum][TVatom_linenum[1]] + m*temp_coordinate_y[tempEQnum][TVatom_linenum[2]];
                  double tempmeta_Z = temp_coordinate_z[tempEQnum][j] + k*temp_coordinate_z[tempEQnum][TVatom_linenum[0]] + l*temp_coordinate_z[tempEQnum][TVatom_linenum[1]] + m*temp_coordinate_z[tempEQnum][TVatom_linenum[2]];
                  fprintf(output, "%-2s\t%17.12lf %17.12lf %17.12lf\n",atomname[tempEQnum][j], tempmeta_X, tempmeta_Y, tempmeta_Z);
                  if( gjfgen_option == 1 ){
                    fprintf(gjf, " %-2s\t%17.12lf %17.12lf %17.12lf\n",atomname[tempEQnum][j], tempmeta_X, tempmeta_Y, tempmeta_Z);
                  }
/*                  if( forDOSinpgen_option == 1 ){ // ( DOS inp #2 BEGIN)
                    fprintf(forDOSinp, "%17.12lf %17.12lf %17.12lf %d\n",tempmeta_X, tempmeta_Y, tempmeta_Z, forDOSinp_Species_num_list[tempnum]);
                  } // ( DOS inp #2 END)
                  tempnum++;*/
                }
              }
            }
          }
        }
        if( gjfgen_option == 1 ){
          fprintf(gjf, " Tv\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[0]] * (xyzframescale_x + 1), temp_coordinate_y[tempEQnum][TVatom_linenum[0]] * (xyzframescale_x + 1), temp_coordinate_z[tempEQnum][TVatom_linenum[0]] * (xyzframescale_x + 1));
          fprintf(gjf, " Tv\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[1]] * (xyzframescale_y + 1), temp_coordinate_y[tempEQnum][TVatom_linenum[1]] * (xyzframescale_y + 1), temp_coordinate_z[tempEQnum][TVatom_linenum[1]] * (xyzframescale_y + 1));
          fprintf(gjf, " Tv\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[2]] * (xyzframescale_z + 1), temp_coordinate_y[tempEQnum][TVatom_linenum[2]] * (xyzframescale_z + 1), temp_coordinate_z[tempEQnum][TVatom_linenum[2]] * (xyzframescale_z + 1));
          fclose(gjf);
        }
        if( supercellgen_option == 1 ){
          sprintf(spfilename, "EQCS_%s_SPcell_EQ%d.gjf", JOBNAME, spEQnum);
          supercellgjf = fopen(spfilename,"w"); //open the output file
          if(supercellxyz == NULL) { //failed to open, return NULL
            printf( "failed to open %s!\n", spfilename);
            return -1;
          }
          fprintf(supercellgjf, "%%chk=%s_EQ%d.chk\n# pbepbe/6-31g/auto\n\n%s_SPcell_EQ%d\n\n0 1\n", JOBNAME, spEQnum, JOBNAME, spEQnum);
          for(j = 0;j < meta_atoms;j++){
            if( j != TVatom_linenum[0] && j != TVatom_linenum[1] && j != TVatom_linenum[2] ){
              fprintf(supercellEQlist, "%-2s\t%17.12lf %17.12lf %17.12lf\n",atomname[tempEQnum][j], temp_coordinate_x[tempEQnum][j], temp_coordinate_y[tempEQnum][j], temp_coordinate_z[tempEQnum][j]);
              fprintf(supercellxyz, "%-2s\t%17.12lf %17.12lf %17.12lf\n",atomname[tempEQnum][j], temp_coordinate_x[tempEQnum][j], temp_coordinate_y[tempEQnum][j], temp_coordinate_z[tempEQnum][j]);
              fprintf(supercellgjf, " %-2s\t%17.12lf %17.12lf %17.12lf\n",atomname[tempEQnum][j], temp_coordinate_x[tempEQnum][j], temp_coordinate_y[tempEQnum][j], temp_coordinate_z[tempEQnum][j]);
            }
          }
          
          for(k = 2;k <= sp_scale_factor;k++){ // added 9722
            for(j = 0;j < meta_atoms;j++){
              if( j != TVatom_linenum[0] && j != TVatom_linenum[1] && j != TVatom_linenum[2] ){
                double tempmeta_X = temp_coordinate_x[tempEQnum][j] + (k-1)*temp_coordinate_x[tempEQnum][TVatom_linenum[shortest_TVnum[tempEQnum]]];
                double tempmeta_Y = temp_coordinate_y[tempEQnum][j] + (k-1)*temp_coordinate_y[tempEQnum][TVatom_linenum[shortest_TVnum[tempEQnum]]];
                double tempmeta_Z = temp_coordinate_z[tempEQnum][j] + (k-1)*temp_coordinate_z[tempEQnum][TVatom_linenum[shortest_TVnum[tempEQnum]]];
                fprintf(supercellEQlist, "%-2s\t%17.12lf %17.12lf %17.12lf\n",atomname[tempEQnum][j], tempmeta_X, tempmeta_Y, tempmeta_Z);
                fprintf(supercellxyz, "%-2s\t%17.12lf %17.12lf %17.12lf\n",atomname[tempEQnum][j], tempmeta_X, tempmeta_Y, tempmeta_Z);
                fprintf(supercellgjf, " %-2s\t%17.12lf %17.12lf %17.12lf\n",atomname[tempEQnum][j], tempmeta_X, tempmeta_Y, tempmeta_Z);
              }
            }
          }
          fprintf(supercellxyz, "\n"); // super cell xyz file and EQ list
          if( shortest_TVnum[tempEQnum] == 0 ){
            fprintf(supercellEQlist, "TV\t%17.12lf %17.12lf %17.12lf\n", sp_scale_factor*temp_coordinate_x[tempEQnum][TVatom_linenum[0]], sp_scale_factor*temp_coordinate_y[tempEQnum][TVatom_linenum[0]], sp_scale_factor*temp_coordinate_z[tempEQnum][TVatom_linenum[0]]);
            fprintf(supercellEQlist, "TV\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[1]], temp_coordinate_y[tempEQnum][TVatom_linenum[1]], temp_coordinate_z[tempEQnum][TVatom_linenum[1]]);
            fprintf(supercellEQlist, "TV\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[2]], temp_coordinate_y[tempEQnum][TVatom_linenum[2]], temp_coordinate_z[tempEQnum][TVatom_linenum[2]]);
            fprintf(supercellgjf, " Tv\t%17.12lf %17.12lf %17.12lf\n", sp_scale_factor*temp_coordinate_x[tempEQnum][TVatom_linenum[0]], sp_scale_factor*temp_coordinate_y[tempEQnum][TVatom_linenum[0]], sp_scale_factor*temp_coordinate_z[tempEQnum][TVatom_linenum[0]]);
            fprintf(supercellgjf, " Tv\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[1]], temp_coordinate_y[tempEQnum][TVatom_linenum[1]], temp_coordinate_z[tempEQnum][TVatom_linenum[1]]);
            fprintf(supercellgjf, " Tv\t%17.12lf %17.12lf %17.12lf\n\n\n", temp_coordinate_x[tempEQnum][TVatom_linenum[2]], temp_coordinate_y[tempEQnum][TVatom_linenum[2]], temp_coordinate_z[tempEQnum][TVatom_linenum[2]]);
          }else if( shortest_TVnum[tempEQnum] == 1 ){
            fprintf(supercellEQlist, "TV\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[0]], temp_coordinate_y[tempEQnum][TVatom_linenum[0]], temp_coordinate_z[tempEQnum][TVatom_linenum[0]]);
            fprintf(supercellEQlist, "TV\t%17.12lf %17.12lf %17.12lf\n", sp_scale_factor*temp_coordinate_x[tempEQnum][TVatom_linenum[1]], sp_scale_factor*temp_coordinate_y[tempEQnum][TVatom_linenum[1]], sp_scale_factor*temp_coordinate_z[tempEQnum][TVatom_linenum[1]]);
            fprintf(supercellEQlist, "TV\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[2]], temp_coordinate_y[tempEQnum][TVatom_linenum[2]], temp_coordinate_z[tempEQnum][TVatom_linenum[2]]);
            fprintf(supercellgjf, " Tv\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[0]], temp_coordinate_y[tempEQnum][TVatom_linenum[0]], temp_coordinate_z[tempEQnum][TVatom_linenum[0]]);
            fprintf(supercellgjf, " Tv\t%17.12lf %17.12lf %17.12lf\n", sp_scale_factor*temp_coordinate_x[tempEQnum][TVatom_linenum[1]], sp_scale_factor*temp_coordinate_y[tempEQnum][TVatom_linenum[1]], sp_scale_factor*temp_coordinate_z[tempEQnum][TVatom_linenum[1]]);
            fprintf(supercellgjf, " Tv\t%17.12lf %17.12lf %17.12lf\n\n\n", temp_coordinate_x[tempEQnum][TVatom_linenum[2]], temp_coordinate_y[tempEQnum][TVatom_linenum[2]], temp_coordinate_z[tempEQnum][TVatom_linenum[2]]);
          }else if( shortest_TVnum[tempEQnum] == 2 ){
            fprintf(supercellEQlist, "TV\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[0]], temp_coordinate_y[tempEQnum][TVatom_linenum[0]], temp_coordinate_z[tempEQnum][TVatom_linenum[0]]);
            fprintf(supercellEQlist, "TV\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[1]], temp_coordinate_y[tempEQnum][TVatom_linenum[1]], temp_coordinate_z[tempEQnum][TVatom_linenum[1]]);
            fprintf(supercellEQlist, "TV\t%17.12lf %17.12lf %17.12lf\n", sp_scale_factor*temp_coordinate_x[tempEQnum][TVatom_linenum[2]], sp_scale_factor*temp_coordinate_y[tempEQnum][TVatom_linenum[2]], sp_scale_factor*temp_coordinate_z[tempEQnum][TVatom_linenum[2]]);
            fprintf(supercellgjf, " Tv\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[0]], temp_coordinate_y[tempEQnum][TVatom_linenum[0]], temp_coordinate_z[tempEQnum][TVatom_linenum[0]]);
            fprintf(supercellgjf, " Tv\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[1]], temp_coordinate_y[tempEQnum][TVatom_linenum[1]], temp_coordinate_z[tempEQnum][TVatom_linenum[1]]);
            fprintf(supercellgjf, " Tv\t%17.12lf %17.12lf %17.12lf\n\n\n", sp_scale_factor*temp_coordinate_x[tempEQnum][TVatom_linenum[2]], sp_scale_factor*temp_coordinate_y[tempEQnum][TVatom_linenum[2]], sp_scale_factor*temp_coordinate_z[tempEQnum][TVatom_linenum[2]]);
          }
          fclose(supercellgjf);
        }        
        if( forDOSinpgen_option == 1 ){ // ( DOS inp #3 BEGIN)
          fprintf(forDOSinp, "%%endblock AtomicCoordinatesAndAtomicSpecies\n\nLatticeConstant             1.0 Ang\n%%block LatticeVectors\n	12.71647794 	0.00000000 	0.00000000 \n	0.00000000 	12.90846708 	0.00000000 \n	0.00000000 	0.00000000 	25.81577617 \n%%endblock LatticeVectors\n\n%%block kgrid_Monkhorst_Pack\n %d  0  0  0.\n 0  %d  0  0.\n 0  0  %d  0.\n%%endblock kgrid_Monkhorst_Pack\n\nPAO.EnergyShift         50 meV\n#%%block PAO.Basis\n#  H_GGA    2      0.0000\n#  0    2 P 1 E     6.0000055      0.2924866\n#  6.0027460        2.0506972\n#  1.00000000000000        1.00000000000000\n#  1    1 P 1 E    39.9677761      0.5000000\n#  2.0285088\n#  1.00000000000000\n#%%endblock PAO.Basis\n\n# Spin\nSpinPolarized               .true.\nFixSpin                     .false.\n#TotalSpin                  1.0\n\n# DFT, Grid, SCF\nXC.functional               GGA         # Exchange-correlation functional type\nXC.authors                  PBE         # Particular parametrization of xc func\nMeshCutoff                  %d. Ry     # Equivalent planewave cutoff for the grid\nMaxSCFIterations            128          # Maximum number of SCF iterations per step\n\n", DOSk1, DOSk2, DOSk3, DOScutoff);
          fprintf(forDOSinp, "DM.MixingWeight             0.1         # New DM amount for next SCF cycle\nDM.Tolerance                1.d-5       # Tolerance in maximum difference between input and output DM\nDM.NumberPulay              20          # Number of SCF steps between pulay mixing\n\nSolutionMethod              Diagon      # OrderN or Diagon\n#OccupationFunction          MP          # The occupation function proposed by Methfessel and Paxton\nElectronicTemperature       300 K     # Temp. for Fermi smearing\n\n# Molecular dynamics and relaxations\nMD.TypeOfRun                CG\nMD.NumCGsteps               0\nMD.MaxForceTol              0.02 eV/Ang # Tolerance in maximum force (default is 0.04 ev/Ang)\nMD.VariableCell             .true.      # Added by maki\n\n# Output options\nWriteCoorInitial\nWriteCoorStep\nWriteForces                 .false.\nWriteKpoints                .false.\nWriteEigenvalues            .true.\nWriteKbands                 .false.\nWriteBands                  .false.\nWriteMullikenPop            1           # Write Mulliken Population Analysis\nWriteCoorXmol               .true.\n\n# Options for saving/reading information\nDM.UseSaveDM                .false.     # Use DM Continuation files\nMD.UseSaveXV                .false.     # Use stored positions and velocities\nMD.UseSaveCG                .false.     # Use stored positions and velocities\nSaveRho                     .false.     # Write valence pseudocharge at the mesh\nSaveDeltaRho                .false.     # Write RHOscf-RHOatm at the mesh\nSaveElectrostaticPotential  .false.     # Write the total elect. pot. at the mesh\nSaveTotalPotential          .false.     # Write the total pot. at the mesh\nWriteSiestaDim              .false.     # Write minimum dim to siesta.h and stop\nWriteDenchar                .false.     # Write information for DENCHAR\n\n\n");
          fprintf(forDOSsh, "#!/bin/csh\n\nset JOB=DOS_%s_EQ%d\n\nset DATA=$PWD\nset WORK=/scr/$USER.${JOB}.$$\n\nsetenv subsiesta siesta\nalias cp 'cp -pf'\n\nif ( -d ${WORK} ) rm -rf ${WORK}\nmkdir -m 700 ${WORK}\n\ncp -f ${JOB}.inp ${WORK}\ncp -f *.psf ${WORK}\ncd ${WORK}\n\nmpirun -np %d ${subsiesta} < ${JOB}.inp > ${JOB}.out\n\ncp -f * ${DATA}\ncd ${DATA}\n\nrm -rf ${WORK}\n\n\n", JOBNAME, tempEQnum, DOSnp);
          fclose(forDOSinp);
          fclose(forDOSsh);
        } // ( DOS inp #3 END)
        if( EQlist2_mode == 0 ){ // if original file is not EQ_list_2 (begin)
            /// make scaled EQ_list (begin)
            for(k = 0;k <= xyzframescale_x;k++){
              for(l = 0;l <= xyzframescale_y;l++){
                for(m = 0;m <= xyzframescale_z;m++){
                  for(j = 0;j < atoms+totalTVnum-1;j++){
                    if( j != TVatom_linenum[0] && j != TVatom_linenum[1] && j != TVatom_linenum[2] && k+l+m != 0 ){
                      double tempmeta_X = temp_coordinate_x[tempEQnum][j] + k*temp_coordinate_x[tempEQnum][TVatom_linenum[0]] + l*temp_coordinate_x[tempEQnum][TVatom_linenum[1]] + m*temp_coordinate_x[tempEQnum][TVatom_linenum[2]];
                      double tempmeta_Y = temp_coordinate_y[tempEQnum][j] + k*temp_coordinate_y[tempEQnum][TVatom_linenum[0]] + l*temp_coordinate_y[tempEQnum][TVatom_linenum[1]] + m*temp_coordinate_y[tempEQnum][TVatom_linenum[2]];
                      double tempmeta_Z = temp_coordinate_z[tempEQnum][j] + k*temp_coordinate_z[tempEQnum][TVatom_linenum[0]] + l*temp_coordinate_z[tempEQnum][TVatom_linenum[1]] + m*temp_coordinate_z[tempEQnum][TVatom_linenum[2]];
                      fprintf(scaledEQlist, "%-2s\t%17.12lf %17.12lf %17.12lf\n",atomname[tempEQnum][j], tempmeta_X, tempmeta_Y, tempmeta_Z);
                    }
                  }
                }
              }
            }
          fprintf(scaledEQlist, "TV\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[0]] * (xyzframescale_x + 1), temp_coordinate_y[tempEQnum][TVatom_linenum[0]] * (xyzframescale_x + 1), temp_coordinate_z[tempEQnum][TVatom_linenum[0]] * (xyzframescale_x + 1));
          fprintf(scaledEQlist, "TV\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[1]] * (xyzframescale_y + 1), temp_coordinate_y[tempEQnum][TVatom_linenum[1]] * (xyzframescale_y + 1), temp_coordinate_z[tempEQnum][TVatom_linenum[1]] * (xyzframescale_y + 1));
          fprintf(scaledEQlist, "TV\t%17.12lf %17.12lf %17.12lf\n", temp_coordinate_x[tempEQnum][TVatom_linenum[2]] * (xyzframescale_z + 1), temp_coordinate_y[tempEQnum][TVatom_linenum[2]] * (xyzframescale_z + 1), temp_coordinate_z[tempEQnum][TVatom_linenum[2]] * (xyzframescale_z + 1));
          if( supercellgen_option == 1 ){
            fprintf(supercellEQlist, "Energy    = %17.12lf ( %17.12lf : %17.12lf)\n", Forced_Energy[tempEQnum] * sp_scale_factor, Bare_Energy[tempEQnum] * sp_scale_factor, Zero_Energy[tempEQnum] * sp_scale_factor);
            for(k = 0;k < totalEQinfonum;k++){
              fprintf(supercellEQlist, "%s",EQinfo[tempEQnum][k]);
            }
          }
          fprintf(scaledEQlist, "Energy    = %17.12lf ( %17.12lf : %17.12lf)\n", Forced_Energy[tempEQnum] * (xyzframescale_x + 1) * (xyzframescale_y + 1) * (xyzframescale_z + 1), Bare_Energy[tempEQnum] * (xyzframescale_x + 1) * (xyzframescale_y + 1) * (xyzframescale_z + 1), Zero_Energy[tempEQnum] * (xyzframescale_x + 1) * (xyzframescale_y + 1) * (xyzframescale_z + 1));

          for(k = 0;k < totalEQinfonum;k++){
            fprintf(scaledEQlist, "%s",EQinfo[tempEQnum][k]);
          }
          SPscaled_nmodenum = ( atoms * sp_scale_factor + totalTVnum ) * 3 - 6;
          scaled_nmodenum = ( atoms * (xyzframescale_x + 1) * (xyzframescale_y + 1) * (xyzframescale_z + 1) + totalTVnum ) * 3 - 6;
          if( EQlist2_mode == 0 ){
            if( supercellgen_option == 1 ){
              fprintf(supercellEQlist, "Normal mode eigenvalues : nmode = %d\n",SPscaled_nmodenum);
            }
            fprintf(scaledEQlist, "Normal mode eigenvalues : nmode = %d\n",scaled_nmodenum);
          }else if( EQlist2_mode == 1 ){
            fprintf(scaledEQlist, "Harmonic frequencies : nmode = %d\n",scaled_nmodenum);
          }
          
          temp_nmodenum = scaled_nmodenum / 5 ;
          for(k = 0;k < temp_nmodenum;k++){
            fprintf(scaledEQlist, "  0.300000000   0.300000000   0.300000000   0.300000000   0.300000000\n");
          }
          SPtemp_nmodenum = SPscaled_nmodenum / 5 ;
          if( supercellgen_option == 1 ){
            for(k = 0;k < SPtemp_nmodenum;k++){
              fprintf(supercellEQlist, "  0.300000000   0.300000000   0.300000000   0.300000000   0.300000000\n");
            }
          }
          temp_nmodenum = scaled_nmodenum % 5 ;
          if( temp_nmodenum > 0){
            fprintf(scaledEQlist, "  0.300000000");
            for(k = 1;k < temp_nmodenum;k++){
              fprintf(scaledEQlist, "   0.300000000");
            }
            fprintf(scaledEQlist, "\n\n");
          }
          SPtemp_nmodenum = SPscaled_nmodenum % 5 ;
          if( SPtemp_nmodenum > 0 && supercellgen_option == 1 ){
            fprintf(supercellEQlist, "  0.300000000");
            for(k = 1;k < SPtemp_nmodenum;k++){
              fprintf(supercellEQlist, "   0.300000000");
            }
            fprintf(supercellEQlist, "\n\n");
          }
        } // if original file is not EQ_list_2 (end)
          cell_Vol[tempEQnum] = temp_coordinate_x[tempEQnum][TVatom_linenum[0]] * temp_coordinate_y[tempEQnum][TVatom_linenum[1]] * temp_coordinate_z[tempEQnum][TVatom_linenum[2]] + temp_coordinate_y[tempEQnum][TVatom_linenum[0]] * temp_coordinate_z[tempEQnum][TVatom_linenum[1]] * temp_coordinate_x[tempEQnum][TVatom_linenum[2]] + temp_coordinate_z[tempEQnum][TVatom_linenum[0]] * temp_coordinate_x[tempEQnum][TVatom_linenum[1]] * temp_coordinate_y[tempEQnum][TVatom_linenum[2]] - temp_coordinate_x[tempEQnum][TVatom_linenum[0]] * temp_coordinate_z[tempEQnum][TVatom_linenum[1]] * temp_coordinate_y[tempEQnum][TVatom_linenum[2]] - temp_coordinate_y[tempEQnum][TVatom_linenum[0]] * temp_coordinate_x[tempEQnum][TVatom_linenum[1]] * temp_coordinate_z[tempEQnum][TVatom_linenum[2]] - temp_coordinate_z[tempEQnum][TVatom_linenum[0]] * temp_coordinate_y[tempEQnum][TVatom_linenum[1]] * temp_coordinate_x[tempEQnum][TVatom_linenum[2]];
          cell_density[tempEQnum] = unit_cell_mass / cell_Vol[tempEQnum] * 1.66054 ; // 1.66054... originates from "atomic mass unit"
          fprintf(cellVDlist, "EQ %*d :     %17.12lf :           %6.3lf\n", totalEQ_digitnum, tempEQnum, cell_Vol[tempEQnum], cell_density[tempEQnum]);
          /// make scaled EQ_list (end)
      } ////////// when crystal or surface mode is on (end)
      printf("\r  runing... (EQ-sorting :EQ %d)",tempEQnum);
      tempEQnum++;
      spEQnum++;
      fprintf(output,"\n");
    }
  }
  if( crystal_surface_mode == 1 ){
    fclose(cellVDlist); // close the scaled_EQ_list file
    if( EQlist2_mode == 1 ){
      fclose(scaledEQlist); // close the scaled_EQ_list file
    }
    if( supercellgen_option == 1 ){
      fclose(supercellEQlist);
    }
    if( cifgen_option == 1 ){
      fclose(ciflist);
    }
  }
  fclose(output); // close the output file
  fclose(fp); // close the original file
  
  // エネルギーの計算とソート
  for(k = 0;k < totalEQnum ;k++){
    sorted_Bare_Energy[k] = Bare_Energy[k];
  }
  for(k = 0;k < totalEQnum ;k++){
    for(m = k+1;m < totalEQnum ;m++){
      if( sorted_Bare_Energy[k] > sorted_Bare_Energy[m] ){
        i = Energy_rank_EQnum[k];
        temp_dnum1 = sorted_Bare_Energy[k];
        Energy_rank_EQnum[k] = Energy_rank_EQnum[m];
        sorted_Bare_Energy[k] = sorted_Bare_Energy[m];
        Energy_rank_EQnum[m] = i;
        sorted_Bare_Energy[m] = temp_dnum1;
      }
    }
  }
  // 9615 moved
  if( energy_sort_mode == 1 && EQlist2_mode == 0 && composition_sorting_option == 0 ){
    if( density_sort_mode == 1 ){
      fprintf(energy_sort,"rank \tEQ# \trel. E(el) (kJ/mol) ( hartree : g/cm^3 )\n");
    }else{
      fprintf(energy_sort,"rank \tEQ# \trel. E(el) (kJ/mol) ( hartree )\n");
    }
    for(k = 0;k < totalEQnum ;k++){
      for(m = 0;m < totalEQnum ;m++){
        if( Energy_rank_EQnum[m] == k ){ // m番目にエネルギーが低いEQ#がkのとき
          if( density_sort_mode == 1 ){
            fprintf(energy_sort," %*d \tEQ %-*d \t%12.7lf ( %17.12lf : %6.3lf )\n", totalEQ_digitnum, k+1, totalEQ_digitnum, Energy_rank_EQnum[k], (sorted_Bare_Energy[k] - sorted_Bare_Energy[0]) * 2625.49962, sorted_Bare_Energy[k], cell_density[Energy_rank_EQnum[k]]); // (9419 added)
          }else{
            fprintf(energy_sort," %*d \tEQ %-*d \t%12.7lf ( %17.12lf )\n", totalEQ_digitnum, k+1, totalEQ_digitnum, Energy_rank_EQnum[k], (sorted_Bare_Energy[k] - sorted_Bare_Energy[0]) * 2625.49962, sorted_Bare_Energy[k]);
          }
          // 23606 added
          if(crystal_surface_mode == 0){
            fprintf(energy_sort_xyz,"%d\n# EQ %d ; %12.7lf ( %17.12lf )\n", atoms, Energy_rank_EQnum[k], (sorted_Bare_Energy[k] - sorted_Bare_Energy[0]) * 2625.49962, sorted_Bare_Energy[k]);
            for(j = 0;j < atoms ;j++){
              fprintf(energy_sort_xyz,"%s    %lf    %lf    %lf\n", atomname[Energy_rank_EQnum[k]][j], temp_coordinate_x[Energy_rank_EQnum[k]][j], temp_coordinate_y[Energy_rank_EQnum[k]][j], temp_coordinate_z[Energy_rank_EQnum[k]][j]);
            }
          }else if(crystal_surface_mode == 1){
            fprintf(energy_sort_xyz,"%d\n# EQ %d ; %12.7lf ( %17.12lf )\n", atoms + atoms * (framescale_x * framescale_y * framescale_z - 1), Energy_rank_EQnum[k], (sorted_Bare_Energy[k] - sorted_Bare_Energy[0]) * 2625.49962, sorted_Bare_Energy[k]);
            for(jj = 0;jj < meta_atoms;jj++){
              if( jj != TVatom_linenum[0] && jj != TVatom_linenum[1] && jj != TVatom_linenum[2]){
                double tempmeta_X = temp_coordinate_x[Energy_rank_EQnum[k]][jj];
                double tempmeta_Y = temp_coordinate_y[Energy_rank_EQnum[k]][jj];
                double tempmeta_Z = temp_coordinate_z[Energy_rank_EQnum[k]][jj];
                fprintf(energy_sort_xyz, "%-2s\t%17.12lf %17.12lf %17.12lf\n",atomname[Energy_rank_EQnum[k]][jj], tempmeta_X, tempmeta_Y, tempmeta_Z);
              }
            }
            for(kk = 0;kk <= framescale_x;kk++){
              for(ll = 0;ll <= framescale_y;ll++){
                for(mm = 0;mm <= framescale_z;mm++){
                  for(jj = 0;jj < meta_atoms;jj++){
                    if( jj != TVatom_linenum[0] && jj != TVatom_linenum[1] && jj != TVatom_linenum[2] && kk+ll+mm != 0 ){
                      double tempmeta_X = temp_coordinate_x[Energy_rank_EQnum[k]][jj] + kk*temp_coordinate_x[Energy_rank_EQnum[k]][TVatom_linenum[0]] + ll*temp_coordinate_x[Energy_rank_EQnum[k]][TVatom_linenum[1]] + mm*temp_coordinate_x[Energy_rank_EQnum[k]][TVatom_linenum[2]];
                      double tempmeta_Y = temp_coordinate_y[Energy_rank_EQnum[k]][jj] + kk*temp_coordinate_y[Energy_rank_EQnum[k]][TVatom_linenum[0]] + ll*temp_coordinate_y[Energy_rank_EQnum[k]][TVatom_linenum[1]] + mm*temp_coordinate_y[Energy_rank_EQnum[k]][TVatom_linenum[2]];
                      double tempmeta_Z = temp_coordinate_z[Energy_rank_EQnum[k]][jj] + kk*temp_coordinate_z[Energy_rank_EQnum[k]][TVatom_linenum[0]] + ll*temp_coordinate_z[Energy_rank_EQnum[k]][TVatom_linenum[1]] + mm*temp_coordinate_z[Energy_rank_EQnum[k]][TVatom_linenum[2]];
                      fprintf(energy_sort_xyz, "%-2s\t%17.12lf %17.12lf %17.12lf\n",atomname[Energy_rank_EQnum[k]][jj], tempmeta_X, tempmeta_Y, tempmeta_Z);
                    }
                  }
                }
              }
            }
          }
          fprintf(energy_sort_xyz,"\n");
          // printf("\n%d", k);
          break;
        }
      }
    }
    fclose(energy_sort);
    fclose(energy_sort_xyz);
  }
  // for density sorting (9419 added)
  if(density_sort_mode == 1 && crystal_surface_mode == 1){
    for(k = 0;k < totalEQnum ;k++){
      sorted_cell_density[k] = cell_density[k];
    }
    for(k = 0;k < totalEQnum ;k++){
      for(m = k+1;m < totalEQnum ;m++){
        if( sorted_cell_density[k] > sorted_cell_density[m] ){
          i = Density_rank_EQnum[k];
          temp_dnum1 = sorted_cell_density[k];
          Density_rank_EQnum[k] = Density_rank_EQnum[m];
          sorted_cell_density[k] = sorted_cell_density[m];
          Density_rank_EQnum[m] = i;
          sorted_cell_density[m] = temp_dnum1;
        }
      }
    }
    fprintf(density_sort,"rank \tEQ No. \td (g/cm^3)\trel. E(el) (kJ/mol)  \tE(el) (hartree)\n");
    for(k = 0;k < totalEQnum ;k++){
      for(m = 0;m < totalEQnum ;m++){
        if( Density_rank_EQnum[m] == k ){
          fprintf(density_sort," %*d \tEQ %-*d \t%6.3lf    \t%12.7lf    \t%17.12lf\n", totalEQ_digitnum, k+1, totalEQ_digitnum, Density_rank_EQnum[k], sorted_cell_density[k], (Bare_Energy[Density_rank_EQnum[k]] - sorted_Bare_Energy[0]) * 2625.49962, Bare_Energy[Density_rank_EQnum[k]]);
          break;
        }
      }
    }
    fclose(density_sort);
  }

  if( gjfgen_option == 1 ){ // 9613 added
    sprintf(command, "mkdir -p EQCS_%s_gjflist; mv EQCS_%s_EQ*.gjf EQCS_%s_gjflist/", JOBNAME, JOBNAME, JOBNAME);
    system(command);
  }
  if( cifgen_option == 1 ){ // 9613 added
    sprintf(command, "mkdir -p EQCS_%s_ciflist; mv EQCS_%s_EQ*.cif EQCS_%s_ciflist/", JOBNAME, JOBNAME, JOBNAME);
    system(command);
  }
  if( supercellgen_option == 1 ){ // 9613 added
    sprintf(command, "mkdir -p EQCS_%s_gjflist; mv EQCS_%s_SPcell_EQ*.gjf EQCS_%s_gjflist/", JOBNAME, JOBNAME, JOBNAME);
    system(command);
  }
  if( forDOSinpgen_option == 1 ){ // 9821 added , 9822 modified
    sprintf(forDOSshfilename, "DOS_%s_exe.sh", JOBNAME);
    forDOSsh = fopen(forDOSshfilename,"w"); //open the output file
    if(forDOSsh == NULL) { //failed to open, return NULL
      printf( "failed to open %s!\n", forDOSshfilename);
      return -1;
    }
    fprintf(forDOSsh, "#!/bin/csh\n\nset JOB=DOS_%s_EQ\nset DATA=$PWD\nset i=1\nwhile( $i <= 218 )\nbsub -n %d -R %sspan[hosts=1]%s -m %sgr11 gr12 gr13%s %s${DATA}/${JOB}${i}.sh%s\n@ i++\nend\n\n\n", JOBNAME, DOSnp, WQ, WQ, WQ, WQ, WQ, WQ);
    fclose(forDOSsh);
    sprintf(command, "mkdir -p DOS_%s_inplist; mv DOS_%s_EQ*.inp DOS_%s_inplist/", JOBNAME, JOBNAME, JOBNAME);
    system(command);
    sprintf(command, "chmod u+x DOS_%s_exe.sh; mv DOS_%s_exe.sh DOS_%s_inplist/", JOBNAME, JOBNAME, JOBNAME);
    system(command);
    sprintf(command, "chmod u+x DOS_%s_EQ*.sh; mv DOS_%s_EQ*.sh DOS_%s_inplist/", JOBNAME, JOBNAME, JOBNAME);
    system(command);
  }
  
  if( totalEQnum > 0 && cytoscapelistgen_option == 1 && composition_sorting_option == 0 ){
    fprintf(cytoscapeEQlist, "EQ,energy(ratio)\n");
    for(j = 0;j < totalEQnum ;j++){
      fprintf(cytoscapeEQlist,"%d,%10.6lf\n", j, (Bare_Energy[j] - sorted_Bare_Energy[0]) / (sorted_Bare_Energy[totalEQnum-1] - sorted_Bare_Energy[0]));
    }
    fclose(cytoscapeEQlist);
  }
  
  i = 0;
  printf("\r   the LOWEST Energy EQ : %d / EQ_MINenergy = %17.12lf\n   the HIGHEST Energy EQ : %d / EQ_MAXenergy = %17.12lf\n", Energy_rank_EQnum[0], sorted_Bare_Energy[0], Energy_rank_EQnum[totalEQnum-1], sorted_Bare_Energy[totalEQnum-1]);
  if(density_sort_mode == 1 && crystal_surface_mode == 1){
    printf("   the LOWEST Density EQ : %d, the HIGHEST Density EQ : %d\n", Density_rank_EQnum[0], Density_rank_EQnum[totalEQnum-1]);
  }
  printf("  *EQ list sorting was completed. ");
  
//////////////////////////////////////// EQlist END 
////////////////////////////////////////################################################################################################################
////////////////////////////////////////################################################################################################################
////////////////////////////////////////################################################################################################################
////////////////////////////////////////################################################################################################################
////////////////////////////////////////################################################################################################################
////////////////////////////////////////################################################################################################################
//////////////////////////////////////// EQlist END 

/////////////////////////////// MAIN SECTION ///////////////////////////////

  if( composition_sorting_option == 1 ){ ///// ##### composition sorting mode is on (begin)
  printf("\n");
  printf("\r  runing... (bond check)");
  i = 1;
  int bonding_checker[totalEQnum][atoms+1][atoms+1]; /// bonding_checker is a bond schecking array !! ( if an element is "1", the two atom (No.j & No.k) have a covalent bond ) (default number of element is "0")
    for(j = 0;j < totalEQnum;j++){ /// intialize "bonding_checker"
      for(k = 0;k < atoms ;k++){
        for(m = 0;m < atoms ;m++){
          bonding_checker[j][k][m] = 0;
        }
      }
    }
  
  int bonding_groups[atoms+1]; /// bonding_checker is a composition sorting array !! ( grouping atoms / separate into molecules ) (default number of element is "0")
    for(k = 0;k < atoms ;k++){ /// intialize "bonding_groups"
      bonding_groups[k] = k;
    }
  int bonding_checking_sum;
  
  char tempatomname1[N],tempatomname2[N];
  char tempcomposition[totalEQnum+1][atoms+1][N], temp_tempcomposition[N];
  double temp_distance[totalEQnum][atoms+1][atoms+1];
  
  char extra[100][L], smiles_list[totalEQnum+1][L], inchi_list[totalEQnum+1][L], babel_version[N]; // (8c26added)
  int smiles_compositionnum[totalEQnum+1], inchi_compositionnum[totalEQnum+1], Iseq_compositionnum[totalEQnum+1], /*babel_compositionnum[totalEQnum+1],*/ smiles_composition_groupnum[totalEQnum+1], inchi_composition_groupnum[totalEQnum+1], total_smiles_compositionnum = 0, total_inchi_compositionnum = 0; // (8c26added)
  int totalcompositionnum = 0, compositionnum[totalEQnum+1], tempcompositionpoint, true_compositionnum[totalEQnum+1];
    for(k = 0;k < totalEQnum ;k++){ /// intialize "bonding_groups"
      compositionnum[k] = k;
      smiles_compositionnum[k] = k;
      inchi_compositionnum[k] = k;
      // babel_compositionnum[k] = k;
    }
  if( no_babel_mode == 0 ){ // Open Babel //
    system("babel -V | grep Babel | cut -d ' ' -f -3 > tempOpenBabel.log");
    fp = fopen("tempOpenBabel.log","r"); // open the original file to read
    fgets(extra[0], N, fp);
    if(extra[0] == "") { // failed to open, return NULL
      printf("Open Babel is not available... (use [-n/N] option)\n");
      return -1;
    }
    fgets(str, N, fp);
    sprintf(babel_version,"%s",str);
    fclose(fp);
    remove("tempOpenBabel.log");
  }// Open Babel //
  
  fp = fopen(original_EQlist_fname,"r"); // open the original file to read  <<2nd time>>
  if(fp == NULL) { // failed to open, return NULL
    printf("failed to open %s!\n", original_EQlist_fname);
    return -1;
  }
  
  sprintf(filename, "tempEQcompositionLOG_%s.log", JOBNAME);
  output = fopen(filename,"w"); // open the output file
  if(output == NULL) { // failed to open, return NULL
    printf("failed to open %s!\n", filename);
    return -1;
  }
  
  sprintf(sublogname, "tempEQcomposition_bond_conections_%s.log", JOBNAME);
  sublog = fopen(sublogname,"w"); // open the output file
  if(sublogname == NULL) { // failed to open, return NULL
    printf("failed to open %s!\n", sublogname);
    return -1;
  }

  ////////////// calculation section <<begin>>
  printf("\r  runing... (distance calculation)");
  for(i = 0;i < totalEQnum;i++){ // start of "for loop" for all EQs
    for(j = 0;j < atoms;j++){ // ## temporary[0] is the line of "# Geometry of EQ X, SYMMETRY = xxx  " /// loop for the main atom //// *** ATTENTION "j < k" !!!!! ***
      for(k = j+1;k < atoms;k++){ /// loop for the comparing atoms
        double tempmeta_X = temp_coordinate_x[i][j] - temp_coordinate_x[i][k];
        double tempmeta_Y = temp_coordinate_y[i][j] - temp_coordinate_y[i][k];
        double tempmeta_Z = temp_coordinate_z[i][j] - temp_coordinate_z[i][k];
        double meta_X = pow(tempmeta_X, 2.0);
        double meta_Y = pow(tempmeta_Y, 2.0);
        double meta_Z = pow(tempmeta_Z, 2.0);
        double temproot = meta_X + meta_Y + meta_Z;
        temp_distance[i][j][k] = sqrt(temproot);
      }
    }
  } // end of "for loop" for all EQs
  if( no_babel_mode == 0 ){
    printf("\r  runing... (generating sequences)");
    if( one_by_one_mode == 0 ){
      ///// SMILES & InChI part <<begin>> // (8c26added) /// #babel
      sprintf(command, "babel -i xyz %s -o can tempEQCS_%s.cansmi > 000_000_nu00", xyzfilename, JOBNAME);
      system(command);
      sprintf(command, "cat tempEQCS_%s.cansmi | awk '{ print $1 }' > EQCS_%s.cansmi", JOBNAME, JOBNAME);
      system(command);
      sprintf(tempfname, "tempEQCS_%s.cansmi", JOBNAME);
      remove(tempfname); // remove the temporary canonical SMILES file
      sprintf(command, "babel -i xyz %s -o inchi EQCS_%s.inchi > 000_000_nu00", xyzfilename, JOBNAME);
      system(command);
    // for SMILES
      sprintf(tempfname, "EQCS_%s.cansmi", JOBNAME);
      temp_inputfile = fopen(tempfname,"r"); // get smiles data (open a temporary input file)
      for(j = 0;j < totalEQnum;j++){
        fgets(extra[0], N, temp_inputfile);
        sscanf(extra[0],"%s\n", smiles_list[j]);
        for(k = 0;k < j;k++){ /// #compare inchi sequences
          if( strcmp( smiles_list[j], smiles_list[k] ) == 0 ){
            smiles_compositionnum[j] = k;
            break;
          }
        }
      }
      fclose(temp_inputfile);
      remove(tempfname); // remove the temporary input file
    // for InChI
      sprintf(tempfname, "EQCS_%s.inchi", JOBNAME);
      temp_inputfile = fopen(tempfname,"r"); // get inchi data (open a temporary input file)
      for(j = 0;j < totalEQnum;j++){
        fgets(extra[0], N, temp_inputfile);
        sscanf(extra[0],"%s\n", inchi_list[j]);
        for(k = 0;k < j;k++){ /// #compare inchi sequences
          if( strcmp( inchi_list[j], inchi_list[k] ) == 0 ){
            inchi_compositionnum[j] = k;
            break;
          }
        }
      }
      fclose(temp_inputfile);
      remove(tempfname); // remove the temporary input file
    }else if( one_by_one_mode == 1 ){ // for POTATO (9417added)
      for(j = 0;j < totalEQnum;j++){
        sprintf(tempfname, "tempEQCS_%s.xyz", JOBNAME);
        temp_babel = fopen(tempfname,"w");
        fprintf(temp_babel, "%d\n", atoms);
        for(k = 0;k < atoms;k++){
          fprintf(temp_babel, "%s %17.12lf %17.12lf %17.12lf", atomname[j][k], temp_coordinate_x[j][k]);
        }
        fclose(temp_babel);
        sprintf(command, "babel -i xyz %s -o can EQ%d_%s.cansmi > 000_000_nu00", tempfname, j, JOBNAME);
        system(command);
        sprintf(command, "babel -i xyz %s -o inchi EQ%d_%s.inchi > 000_000_nu00", tempfname, j, JOBNAME);
        system(command);
        remove(tempfname);
        sprintf(tempfname, "EQ%d_%s.cansmi", j, JOBNAME);
        temp_babel = fopen(tempfname,"r"); // get smiles data (open a temporary input file)
        fgets(extra[0], N, temp_babel); // "Are there cases that the sequence length is longer than N?"
        sscanf(extra[0],"%s # %s\n",smiles_list[j],extra[1]);
        fclose(temp_babel);
        remove(tempfname); // remove the temporary input file
        sprintf(tempfname, "EQ%d_%s.inchi", j, JOBNAME);
        temp_babel = fopen(tempfname,"r"); // get inchi data (open a temporary input file)
        fgets(extra[0], N, temp_babel);
        sscanf(extra[0],"%s\n",inchi_list[j]);
        fclose(temp_babel);
        remove(tempfname); // remove the temporary input file
        for(k = 0;k < j;k++){ /// #compare smiles sequences
  //          if( strcount( smiles_list[j],smiles_list[k] ) == 1 && strcount( smiles_list[k],smiles_list[j] ) == 1 ){
          if( strcmp( inchi_list[j],inchi_list[k] ) == 0 ){
            smiles_compositionnum[j] = k;
            break;
          }
        }
        for(k = 0;k < j;k++){ /// #compare inchi sequences
  //          if( strcount( inchi_list[j],inchi_list[k] ) == 1 && strcount( inchi_list[k],inchi_list[j] ) == 1 ){
          if( strcmp( inchi_list[j],inchi_list[k] ) == 0 ){
            inchi_compositionnum[j] = k;
            break;
          }
        }
      }
    }
  }
      ///// SMILES & InChI part <<end>> // (8c26added)
  tempEQnum = 0;
  tempnum = 1;
  ////////////// calculation section <<end>>
  
  ////////////// grouping section for SMILES & InChI <<begin>>
  printf("\r  runing... (sequence grouping)     ");
      i = 0;
      for(j = 0;j < totalEQnum;j++){
        for(k = 0;k < totalEQnum;k++){
          if( smiles_compositionnum[k] == j ){
            smiles_composition_groupnum[k] = i;
            for(m = k+1;m < totalEQnum;m++){
              if( smiles_compositionnum[m] == j ){
                smiles_composition_groupnum[m] = i;
              }
            }
            i++;
            break;
          }
        }
      }
      for(j = 0;j < totalEQnum;j++){
        for(k = 0;k < totalEQnum;k++){
          if( smiles_composition_groupnum[k] == j ){
            total_smiles_compositionnum++;
            break;
          }
        }
      }
      i = 0;
      
      for(j = 0;j < totalEQnum;j++){
        for(k = 0;k < totalEQnum;k++){
          if( inchi_compositionnum[k] == j ){
            inchi_composition_groupnum[k] = i;
            for(m = k+1;m < totalEQnum;m++){
              if( inchi_compositionnum[m] == j ){
                inchi_composition_groupnum[m] = i;
              }
            }
            i++;
            break;
          }
        }
      }
      for(j = 0;j < totalEQnum;j++){
        for(k = 0;k < totalEQnum;k++){
          if( inchi_composition_groupnum[k] == j ){
            total_inchi_compositionnum++;
            break;
          }
        }
      }
      i = 0;
  ////////////// grouping section for SMILES & InChI <<end>>
  
  ////////////////////////////////// script section <<begin>> //////////////////////////////////
  fprintf(output,"*************************************************************************\n   EQ composition sorting LOG (version %s)\n      \"%s\"\n*************************************************************************\n                                            produced by H.N. (%s)\n   the number of EQ      : %d\n   the number of atom    : %d\n   the number of element : %d\n\n=========================================================================\n\n", software_version, JOBNAME, developed_date, totalEQnum, atoms, totalelementnamenum);
  
  fprintf(sublog,"*************************************************************************\n   EQ bond conections LOG (version %s)\n      \"%s\"\n*************************************************************************\n                                            produced by H.N. (%s)\n   the number of EQ      : %d\n   the number of atom    : %d\n   the number of element : %d\n\n=========================================================================\n\n", software_version, JOBNAME, developed_date, totalEQnum, atoms, totalelementnamenum);
  printf("\r  runing... (type grouping)     ");
  while(tempEQnum < totalEQnum){
    fprintf(output, "# EQ %d \n",tempEQnum);

  ///// grouping section <<begin>>
    for(j = 0;j < atoms;j++){ //// *** ATTENTION "j < k" !!!!! *** /// bonding_checker[tempEQnum][j][k] is only valid ( bonding_checker[k][j] = 0 !!! )
      sprintf(atomname2[j][0],"%s",atomname_TVexcluded[0]);
      for(k = j+1;k < atoms;k++){
        sprintf(tempatomname1,"%sX",atomname_TVexcluded[j]);
        sprintf(tempatomname2,"%sX",atomname_TVexcluded[k]);
        double r_sum = radius(tempatomname1) + radius(tempatomname2);
        double r_covalent = 1.25 * r_sum;
        if( temp_distance[tempEQnum][j][k] < r_covalent ){
          bonding_checker[tempEQnum][j][k] = 1; /// the two atom (No.j & No.k) have a covalent bond
          bonding_groups[k] = j; /// input *original* bonding groups
        }
      }
    }
    for(j = 0;j < atoms;j++){ /// redo molecule grouping for listed atoms
      for(k = j+1;k < atoms;k++){
        if( bonding_checker[tempEQnum][j][k] == 1 ){
          bonding_groups[k] = bonding_groups[j]; /// refresh the group list
        }
      }
    }
    for(j = 0;j < atoms;j++){ /// redo molecule grouping for listed atoms
      for(k = j+1;k < atoms;k++){
        for(m = 0;m < atoms;m++){
          if( bonding_checker[tempEQnum][j][m] == 1 && bonding_groups[k] == bonding_groups[m] ){ /// (*)
            bonding_groups[k] = bonding_groups[j];
          }
          if( bonding_checker[tempEQnum][j][k] == 1 && bonding_groups[k] == bonding_groups[m] ){ /// (**)
            bonding_groups[m] = bonding_groups[j];
          }
          if( bonding_checker[tempEQnum][j][m] == 1 && bonding_groups[k] == bonding_groups[m] ){ /// the same casework as (*)
            bonding_groups[k] = bonding_groups[j];
          }
          if( bonding_checker[tempEQnum][j][k] == 1 && bonding_groups[k] == bonding_groups[m] ){ /// the same casework as (**)
            bonding_groups[m] = bonding_groups[j];
          }
          if( bonding_checker[tempEQnum][k][m] == 1 && bonding_groups[j] == bonding_groups[m] ){ /// !! 0~m~atoms % j+1~k~atoms (when k < m)
            bonding_groups[k] = bonding_groups[j];
          }
          if( bonding_checker[tempEQnum][m][k] == 1 && bonding_groups[j] == bonding_groups[m] ){ /// !! 0~m~atoms % j+1~k~atoms (when k > m)
            bonding_groups[k] = bonding_groups[j];
          }
        }
      }
    }
    for(k = 0;k < atoms;k++){
      if( bonding_checker[tempEQnum][k][atoms-1] == 1 ){ /// !! in this loop "k<j"
        bonding_groups[atoms-1] = bonding_groups[k];
        break; // stop by a certain "k"
      }
    }
    for(j = 0;j < tempEQnum;j++){ /// check composition type
      tempnum = 0; tempcompositionpoint = 0;
      for(k = 0;k < atoms;k++){
        for(m = k+1;m < atoms;m++){
          if( bonding_checker[tempEQnum][k][m] == bonding_checker[j][k][m] ){
            tempcompositionpoint++;
          }
          tempnum++;
        }
      }
      if( tempcompositionpoint == tempnum ){ /// if tempcompositionpoint is equal to tempnum, composition number is same
        compositionnum[tempEQnum] = compositionnum[j];
        break;
      }
    }
    if( compositionnum[tempEQnum] == tempEQnum ){
      compositionnum[tempEQnum] = totalcompositionnum;
      totalcompositionnum++;
    }
    tempnum = 1;
    
    
    if ( no_babel_mode == 1 ){
      inchi_composition_groupnum[tempEQnum] = compositionnum[tempEQnum];
    }
    ///// grouping section <<end>> ///// 


    ///// output section <<begin>> ///// 
      if( no_babel_mode == 0 ){
        fprintf(output,"composition type : %d\nsequence type    : %d\nSMILES (Sseq %*d)\n : %s\nInChI  (Iseq %*d)\n : %s\n-------------------------------------------------------------------------\n",compositionnum[tempEQnum]+1,inchi_composition_groupnum[tempEQnum]+1,totalEQ_digitnum,smiles_composition_groupnum[tempEQnum]+1,smiles_list[tempEQnum],totalEQ_digitnum,inchi_composition_groupnum[tempEQnum]+1,inchi_list[tempEQnum]);
      }else{
        fprintf(output,"composition type : %d\n-------------------------------------------------------------------------\n",compositionnum[tempEQnum]+1);
      }
    
    true_compositionnum[tempEQnum] = compositionnum[tempEQnum] + 1;

    i = 1; /// initialize "i"
    for(j = 0;j < atoms;j++){
      for(k = 0;k < atoms;k++){
        if( bonding_groups[k] == j ){ // if No.k atom belongs to No.j bonding group
          fprintf(output,"bonding group : %d\n",i); // (initial value of "i" is 1)
          fprintf(output,"ordered atoms : %s(%d)",atomname[0][k],k+1);
          sprintf(tempcomposition[tempEQnum][i], "%sX", atomname[0][k]); /// add atoms of this No.i bonding group
          for(m = k+1;m < atoms;m++){
            if( bonding_groups[m] == j ){ // if No.m atom belongs to No.j bonding group
              fprintf(output,"%s(%d)", atomname[0][m],m+1);
              sprintf(temp_tempcomposition, "%s", tempcomposition[tempEQnum][i]);
              sprintf(tempcomposition[tempEQnum][i], "%s%sX", temp_tempcomposition, atomname[0][m]); /// add atoms of this No.i bonding group
              tempnum++;
            }
          }
          fprintf(output,"\ncomposition   : ");
          for(m = 0;m < totalelementnamenum;m++){ /// output the compositions for each molecules
            each_atomnum = strcount( tempcomposition[tempEQnum][i], totalelementname[m] ); /// using "strcount" function (original function)
            elementnamelength = strlen(totalelementname[m]);
            sprintf(tempelementname_for_count,"%.*s",elementnamelength-1,totalelementname[m]);
            if( each_atomnum > 1 ){
              fprintf(output,"%s%d",tempelementname_for_count,each_atomnum);
            }else if( each_atomnum == 1 ){
              fprintf(output,"%s",tempelementname_for_count);
            }
          }
          fprintf(output," (%d atoms)\n",tempnum);
          tempnum = 1;
          fprintf(output,"-------------------------------------------------------------------------\n");
          i++;
          break;
        }
      }
    }
    i = 1; /// initialize "i"
    
      /// output into EQcomposition_bond_conections log file (begin)
      if( no_babel_mode == 0 ){
        fprintf(sublog, " # EQ %d \nSMILES : %s\nInChI  : %s\n-------------------------------------------------------------------------\n",tempEQnum,smiles_list[tempEQnum],inchi_list[tempEQnum]);
      }else{
        fprintf(sublog, " # EQ %d \n-------------------------------------------------------------------------\n",tempEQnum);  
      }
        for(j = 0;j < atoms;j++){
          fprintf(sublog,"bonding group (molecule group)         \t: %d\n",bonding_groups[j]);
          fprintf(sublog,"elder number atoms connected to %s(%*d)\t: ",atomname[0][j],totalatoms_digitnum,j+1);
          for(k = j+1;k < atoms;k++){  /// bonding atoms
            if( bonding_checker[tempEQnum][j][k] == 1 ){
              fprintf(sublog,"%s(%d)", atomname2[j][k],k+1);
            }
          }
          fprintf(sublog,"\n");
          for(k = 0;k < atoms;k++){
            if( bonding_checker[tempEQnum][j][k] == 1 ){
              fprintf(sublog, "\t%s(%*d) - %s(%*d)\t= %13.12f\n", atomname[0][j],totalatoms_digitnum, j+1, atomname2[j][k],totalatoms_digitnum, k+1, temp_distance[tempEQnum][j][k]);
            }
          }
          fprintf(sublog,"-------------------------------------------------------------------------\n");
        }
      /// output into EQcomposition_bond_conections log file (end)
      
    ///// output section <<end>>
    if( EQlist2_mode == 0 ){
      fprintf(output,"%s\n=========================================================================\n\n",EQinfo[tempEQnum][0]);
      fprintf(sublog,"%s\n=========================================================================\n\n",EQinfo[tempEQnum][0]);
    }else if( EQlist2_mode == 1 ){
      for(k = 0;k < 3;k++){
        fprintf(output,"%s",EQinfo[tempEQnum][k]);
        fprintf(sublog,"%s",EQinfo[tempEQnum][k]);
      }
      fprintf(output,"%s\n=========================================================================\n\n",EQinfo[tempEQnum][3]);
      fprintf(sublog,"%s\n=========================================================================\n\n",EQinfo[tempEQnum][3]);
    }
    tempEQnum++; /// to NEXT EQs ...
    for(j = 0;j < atoms;j++){ /// redo molecule grouping for listed atoms
      bonding_groups[j] = j; /// initialize the group list
    }
  }
  fprintf(output,"%s  Values of covalent radius are taken from \n\t\"B. Cordero, et al., Dalton Trans. 2008, 0, 2832-2838\"\n  Thresholds are determined by values of \"1.25*(R_i+R_j)\"\n  (Some values of covalent radius are modified; H,Si)\n",asterisk);
  fprintf(sublog,"%s  Values of covalent radius are taken from \n\t\"B. Cordero, et al., Dalton Trans. 2008, 0, 2832-2838\"\n  Thresholds are determined by values of \"1.25*(R_i+R_j)\"\n  (Some values of covalent radius are modified; H,Si)\n",asterisk);
  if( no_babel_mode == 0 ){
    fprintf(output,"  SMILES and InChI sequences are generated by %s",babel_version);
    fprintf(sublog,"  SMILES and InChI sequences are generated by %s",babel_version);
  }
  fprintf(output,"%s",asterisk);
  fprintf(sublog,"%s",asterisk);
  fclose(output); // close the temp output file
  fclose(sublog); // close the temp output (sublog) file
  fclose(fp); // close the original file
  ///////////////////////////////////// script section <<end>> ////////////////////////////////////

  ////////////////////////////////// re-output section <<bigin>> //////////////////////////////////
  printf("\r  runing... (output)");
  char tempLOGfilename[N], composition_str[N], tempcomposition_list[totalcompositionnum][atoms][N], tempatomname[N];
  int tempatomnamelength;
  int bonding_groups_num[totalcompositionnum];
  int totalcomposition_digitnum = digit(totalcompositionnum);
  int temp_babel_sequence_num_list[totalEQnum],temp_babel_sequence_num_listnum; // added 8c31
  int line_num = 1;
  
  sprintf(tempLOGfilename, "tempEQcompositionLOG_%s.log", JOBNAME);
  fp = fopen(tempLOGfilename,"r"); // open the output file
  if(fp == NULL) { // failed to open, return NULL
    printf("failed to open %s!\n", tempLOGfilename);
    return -1;
  }
  j = 0; // "j" is composition type number
  k = 0; // "k" is the number of bonding group of "composition type j"
  while(fgets(str, N, fp) != NULL){ // make tempcompositionlist
    sprintf(composition_str, "composition type : %d\n",j+1); ///*** This string is variable !!!  ***/// <<< CAUTION >>>
    if( strstr( str , composition_str ) != NULL ){ // if "composition type" found in "tempEQcompositionLOG_%.*s.log" file
      while( strstr( str , "======" ) == NULL ){ // until "=" appears
        fgets(str, N, fp);
        if( strstr( str , "composition" ) != NULL ){ // if "composition" appears
          sscanf(str,"composition   : %s (%d atoms)\n",tempcomposition_list[j][k],&tempnum);
          k++;
        }
      }
      bonding_groups_num[j] = k; // keep the number of bonding groups for each composition type
      k = 0; // initialize "k"
      j++;
    }
  }
  tempnum = 0;
  fclose(fp);
  
  sprintf(reinput_filename, "tempEQcompositionLOG_%s.log", JOBNAME);
  fp = fopen(reinput_filename,"r"); // open the output file
  if(fp == NULL) { // failed to open, return NULL
    printf("failed to open %s!\n", reinput_filename);
    return -1;
  }
  sprintf(reoutput_filename, "EQCS_%s_cLOG.log", JOBNAME);
  re_output = fopen(reoutput_filename,"w"); // open the output file
  if(re_output == NULL) { // failed to open, return NULL
    printf("failed to open %s!\n", reoutput_filename);
    return -1;
  }

  i = 0; /// initialize "i" as "0" !!! (9b27 modified)
  while(fgets(str, N, fp) != NULL){
    fprintf(re_output,"%s",str);
    if( line_num == 8 ){
      fprintf(re_output,"                         :");
      for(k = 0;k < totalelementnamenum ;k++){
        tempatomnamelength = strlen(totalelementname[k]);
        sprintf(tempatomname,"%.*s",tempatomnamelength-1,totalelementname[k]);
        fprintf(re_output," %s", tempatomname);
      }
      if( no_babel_mode == 0 ){
        fprintf(re_output,"\n\n   the number of SMILES sequence : %d",total_smiles_compositionnum);
        fprintf(re_output,"\n   the number of InChI sequence  : %d",total_inchi_compositionnum);
      }
      fprintf(re_output,"\n   the number of composition     : %d\n\n", totalcompositionnum);
      fprintf(re_output,"-------------------------------------------------------------------------\n\n   << list of compositions >>\n\n");
      for(j = 0;j < totalcompositionnum;j++){ // list up numbers of EQ whose composition type is "j"
        fprintf(re_output,"   type %d EQ\t:",j+1);
        for(k = 0;k < totalEQnum;k++){
          if( compositionnum[k] == j ){
            fprintf(re_output," %d",k);
            i++;
          }
        }
        fprintf(re_output," (%d EQ(s))\n             \t:",i);
        i = 0;
        tempnum = 0;
        for(k = 0;k < bonding_groups_num[j];k++){
          if( tempcomposition_list[j][k] == NULL ){
            break;
          }else{
            fprintf(re_output," %s", tempcomposition_list[j][k]);
            tempnum++;
          }
        }
        fprintf(re_output," (%d molecule(s))\n",tempnum);
      }
      tempnum = 0;
      fprintf(re_output,"\n"); // added 8c31
      fprintf(re_output,"-------------------------------------------------------------------------\n\n");
      if( no_babel_mode == 0 ){
        fprintf(re_output,"   << list of sequences >>\n\n");
        for(j = 0;j < total_inchi_compositionnum;j++){ // list up numbers of EQ whose sequence type is "j" // added 8c31
          fprintf(re_output,"   Iseq %*d EQ\t:",totalcomposition_digitnum,j+1);
          for(k = 0;k < totalEQnum;k++){
            if( inchi_composition_groupnum[k] == j ){
              fprintf(re_output," %d",k);
              i++;
            }
          }
          fprintf(re_output," (%d EQ(s))\n",i);
          i = 0;
        }
      }
    }
    line_num++;
  }
  i = 1; /// initialize "i" // added 8c31 (option?)
  fclose(re_output); // close the output file
  fclose(fp);
  remove(reinput_filename);
  line_num = 1; /// initialize "i"

  sprintf(reinput_filename, "tempEQcomposition_bond_conections_%s.log", JOBNAME);
  fp = fopen(reinput_filename,"r"); // open the output file
  if(fp == NULL) { // failed to open, return NULL
    printf("failed to open %s!\n", reinput_filename);
    return -1;
  }
  sprintf(reoutput_filename, "EQCS_%s_bLOG.log", JOBNAME);
  re_sublog = fopen(reoutput_filename,"w"); // open the output file
  if(re_sublog == NULL) { // failed to open, return NULL
    printf("failed to open %s!\n", reoutput_filename);
    return -1;
  }

  while(fgets(str, N, fp) != NULL){
    fprintf(re_sublog,"%s",str);
    if( line_num == 8 ){
      fprintf(re_sublog,"                         :");
      for(k = 0;k < totalelementnamenum ;k++){
        tempatomnamelength = strlen(totalelementname[k]);
        sprintf(tempatomname,"%.*s",tempatomnamelength-1,totalelementname[k]);
        fprintf(re_sublog," %s", tempatomname);
      }
      if( no_babel_mode == 0 ){
        fprintf(re_output,"\n\n   the number of SMILES sequence : %d",total_smiles_compositionnum);
        fprintf(re_output,"\n   the number of InChI sequence  : %d",total_inchi_compositionnum);
      }
      fprintf(re_sublog,"\n   the number of composition     : %d\n\n", totalcompositionnum);
      fprintf(re_sublog,"-------------------------------------------------------------------------\n\n   << list of compositions >>\n\n");
      for(j = 0;j < totalcompositionnum;j++){ // list up numbers of EQ whose composition type is "j"
        fprintf(re_sublog,"   type %d EQ\t:",j+1);
        for(k = 0;k < totalEQnum;k++){
          if( compositionnum[k] == j ){
            fprintf(re_sublog," %d",k);
            i++;
          }
        }
        fprintf(re_sublog," (%d EQ(s))\n             \t:",i);
        i = 0;
        tempnum=0;
        for(k = 0;k < bonding_groups_num[j];k++){
          if( tempcomposition_list[j][k] == NULL ){
            break;
          }else{
            fprintf(re_sublog," %s", tempcomposition_list[j][k]);
            tempnum++;
          }
        }
        fprintf(re_output," (%d molecule(s))\n",tempnum);
      }
      tempnum=0;
      fprintf(re_sublog,"\n"); // added 8c31
      fprintf(re_sublog,"-------------------------------------------------------------------------\n\n");
      if( no_babel_mode == 0 ){
        fprintf(re_sublog,"   << list of sequences >>\n\n");
        for(j = 0;j < total_inchi_compositionnum;j++){ // list up numbers of EQ whose sequence type is "j" // added 8c31
          fprintf(re_sublog,"   Iseq %*d EQ\t:",totalcomposition_digitnum,j+1);
          for(k = 0;k < totalEQnum;k++){
            if( inchi_composition_groupnum[k] == j ){
              fprintf(re_sublog," %d",k);
              i++;
            }
          }
          fprintf(re_sublog," (%d EQ(s))\n",i);
          i = 0;
        }
      }
    }
    line_num++;
  }
  i = 1; /// initialize "i" // added 8c31 (option?)
  line_num = 1; /// initialize "line_num"
  fclose(re_sublog); // close the output file
  fclose(fp);
  remove(reinput_filename);

  //////////// 9615 added (EQ list for Cytoscape) (begin)
    if( cytoscapelistgen_option == 1 && composition_sorting_option == 1 && no_babel_mode == 1 ){
      fprintf(cytoscapeEQlist, "EQ,energy(ratio),type\n");
      for(j = 0;j < totalEQnum ;j++){
        fprintf(cytoscapeEQlist,"%d,%10.6lf,%d\n", j, (Bare_Energy[j] - sorted_Bare_Energy[0]) / (sorted_Bare_Energy[totalEQnum-1] - sorted_Bare_Energy[0]), compositionnum[tempEQnum]+1 );
      }
    }else if( cytoscapelistgen_option == 1 && composition_sorting_option == 1 && no_babel_mode == 0 ){
      fprintf(cytoscapeEQlist, "EQ,energy(ratio),type,Iseq\n");
      for(j = 0;j < totalEQnum ;j++){
        fprintf(cytoscapeEQlist,"%d,%10.6lf,%d,%d\n", j, (Bare_Energy[j] - sorted_Bare_Energy[0]) / (sorted_Bare_Energy[totalEQnum-1] - sorted_Bare_Energy[0]), compositionnum[j]+1, inchi_composition_groupnum[j]+1 );
      }
    }
  //////////// 9615 added (EQ list for Cytoscape) (end)

  //////////////////////////////////// re-output section <<end>> ////////////////////////////////////
  // printf("ok.p2,pg\n",no_babel_mode); //################################################################################################################################################################################################
  
  ///////////////////////////////// TS/PT-sorting section <<begin>> /////////////////////////////////
  char TS_fname[N],PT_fname[N],TS_output_fname[N],PT_output_fname[N],tempTS_output_fname[N],tempPT_output_fname[N];
  int TSTS_on = 0, PTTS_on = 0, totalTSTS_digitnum, totalPTTS_digitnum;
  int tempTSTSnum = 0, tempPTTSnum = 0, totalTSTSnum = 0, totalPTTSnum = 0;
  char PTTS_CONNECTION_EQa[totalEQnum][5],PTTS_CONNECTION_EQb[totalEQnum][5];
  char temp_char[N][N];

  if(EQlist2_mode == 1) {
    sprintf( TS_fname,"%s_TS_list_2.log",JOBNAME);
    fp = fopen( TS_fname,"r"); // open the original file to read
    if(fp == NULL) { // failed to open, return NULL
      printf("%s_TS_list_2.log doesn't exist...\n", JOBNAME);
      sprintf(TS_fname,"%s_TS_list.log",JOBNAME);
    }else{
      TSlist2_mode = 1;
    }
  }else{
    sprintf(TS_fname,"%s_TS_list.log",JOBNAME);
    fp = fopen(TS_fname,"r");// open TS_list to read
  }
  if(fp != NULL) { // if TS_list exists <<*begin*>>
    TSTS_on = 1;
    while(fgets(str, N, fp) != NULL){
      if( strstr( str , "CONNECTION" ) != NULL ){
        tempTSTSnum++;
      }
    }
    totalTSTSnum = tempTSTSnum;
    tempTSTSnum = 0; // initialize "tempTSTSnum"
    fclose(fp);
    
    totalTSTS_digitnum = digit(totalTSTSnum);

    if( totalTSTSnum > 0 ){
      printf("\r  runing... (TS-sorting)");
      char TSTS_CONNECTION_EQa[totalTSTSnum][N], TSTS_CONNECTION_EQb[totalTSTSnum][N]/*, TSinfo[totalTSTSnum][totalEQinfonum][N]*/;
      double TSTS_Forced_Energy[totalTSTSnum], TSTS_Bare_Energy[totalTSTSnum], TSTS_Zero_Energy[totalTSTSnum];
      double TS_Electron_Energy[totalTSTSnum], TS_zero_Gibbs_Energy[totalTSTSnum], TS_RT_Gibbs_Energy[totalTSTSnum]; // for TS_list_2 (9507 added) 
      int TSTS_CONNECTION_EQa_num[totalTSTSnum],TSTS_CONNECTION_EQb_num[totalTSTSnum];
      int TSTS_composition_EQa_num[totalTSTSnum],TSTS_composition_EQb_num[totalTSTSnum];
      int TSTS_DC_num[totalTSTSnum]; // 0 (not to DC) or 1 (connected to DC)
      int TSTS_same_composition_list[totalTSTSnum]; // 0 (not same) or 1 (connected to the same composition)
      for(j = 0;j < totalTSTSnum;j++){
        TSTS_same_composition_list[j] = 0;
        TSTS_DC_num[j] = 0;
      }
      sprintf(tempTS_output_fname, "tempEQcomposition_TSlist_%s.log", JOBNAME);
      output = fopen(tempTS_output_fname,"w"); // open the output file
      if(output == NULL) { // failed to open, return NULL
        printf("failed to open %s!\n", tempTS_output_fname);
        return -1;
      }
      if( no_babel_mode == 0 ){
        fprintf(output,"*************************************************************************\n   EQ composition TS list (version %s)\n      \"%s\"\n*************************************************************************\n                                            produced by H.N. (%s)\n   the number of EQ            : %d\n   the number of atom          : %d\n   the number of element       : %d\n   the number of composition   : %d\n   the number of sequence type : %d\n\n-------------------------------------------------------------------------\n\n   << list of compositions >>\n\n", software_version, JOBNAME, developed_date, totalEQnum, atoms, totalelementnamenum, totalcompositionnum, total_inchi_compositionnum);
      }else{
        fprintf(output,"*************************************************************************\n   EQ composition TS list (version %s)\n      \"%s\"\n*************************************************************************\n                                            produced by H.N. (%s)\n   the number of EQ            : %d\n   the number of atom          : %d\n   the number of element       : %d\n   the number of composition   : %d\n\n\n-------------------------------------------------------------------------\n\n   << list of compositions >>\n\n", software_version, JOBNAME, developed_date, totalEQnum, atoms, totalelementnamenum, totalcompositionnum);
      }
      for(j = 0;j < totalcompositionnum;j++){ // list up numbers of EQ whose composition type is "j"
        fprintf(output,"   type %d EQ\t:",j+1);
        i = 0;
        for(k = 0;k < totalEQnum;k++){
          if( compositionnum[k] == j ){
            fprintf(output," %d",k);
            i++;
          }
        }
        fprintf(output," (%d EQ(s))\n             \t:",i);
        i = 0;
        tempnum=0;
        for(k = 0;k < bonding_groups_num[j];k++){
          if( tempcomposition_list[j][k] == NULL ){
            break;
          }else{
            fprintf(output," %s", tempcomposition_list[j][k]);
            tempnum++;
          }
        }
        fprintf(output," (%d molecule(s))\n",tempnum);
      }
      tempnum=0;
      fprintf(output,"\n"); // added 8c31
      if( no_babel_mode == 0 ){
        fprintf(output,"-------------------------------------------------------------------------\n\n   << list of sequences >>\n\n");
        for(j = 0;j < total_inchi_compositionnum;j++){ // list up numbers of EQ whose sequence type is "j" // added 8c31
          fprintf(output,"   Iseq %*d EQ\t:",totalcomposition_digitnum,j+1);
          for(k = 0;k < totalEQnum;k++){
            if( inchi_composition_groupnum[k] == j ){
              fprintf(output," %d",k);
              i++;
            }
          }
          fprintf(output," (%d EQ(s))\n",i);
          i = 0;
        }
      }
      fprintf(output,"\n=========================================================================\n\n");
      fp = fopen(TS_fname,"r");
      i = 0;
      while(fgets(str, N, fp) != NULL){
        if( strstr( str , "Energy" ) != NULL ){
          if( TSlist2_mode == 0 ){
            if( strstr( str , ":" ) != NULL ){
              sscanf(str, "Energy    = %lf (%lf : %lf)\n", &TSTS_Forced_Energy[tempTSTSnum], &TSTS_Bare_Energy[tempTSTSnum], &TSTS_Zero_Energy[tempTSTSnum]); // scan TS barrier
              if( TSTS_Bare_Energy[tempTSTSnum] > -0.1 ){
                TSTS_Bare_Energy[tempTSTSnum] = TSTS_Forced_Energy[tempTSTSnum];
              }
            }else{
              sscanf(str, "Energy    = %lf\n", &TSTS_Bare_Energy[tempTSTSnum]); // scan TS barrier
            }
          }else if( TSlist2_mode == 1 ){ // added 9507
            sscanf(str, "Relat.Energy             = %lf kJ/mol\n", &TS_Electron_Energy[tempTSTSnum]);
            fgets(str, N, fp);
            sscanf(str, "Relat.Energy (   0.00 K) = %lf kJ/mol\n", &TS_zero_Gibbs_Energy[tempTSTSnum]);
            fgets(str, N, fp);
            sscanf(str, "Relat.Energy ( 298.15 K) = %lf kJ/mol\n", &TS_RT_Gibbs_Energy[tempTSTSnum]);
//            fprintf(output, "# EQ %d (Old EQ %d) Rel.E=%.1f; Rel.G(0K)=%.1f; Rel.G(RT)=%.1f; (%.3f%%)\n", tempEQnum, oldEQnum[tempEQnum], Electron_Energy[tempEQnum], zero_Gibbs_Energy[tempEQnum], RT_Gibbs_Energy[tempEQnum], Boltzmann_Distribution[tempEQnum]);
          }
          if( TSTS_MAXenergy < TSTS_Bare_Energy[tempTSTSnum] ){
            TSTS_MAXenergy = TSTS_Bare_Energy[tempTSTSnum];
          }
        }
        if( strstr( str , "CONNECTION" ) != NULL ){
          sscanf(str, "%s %s %s %s %s\n",temp_char[0],temp_char[1],TSTS_CONNECTION_EQa[tempTSTSnum],temp_char[2],TSTS_CONNECTION_EQb[tempTSTSnum]); // scan "CONNECTION : 0 - ??"
          if( strstr(TSTS_CONNECTION_EQa[tempTSTSnum],"??") != NULL || strstr(TSTS_CONNECTION_EQa[tempTSTSnum],"DC") != NULL || strstr(TSTS_CONNECTION_EQb[tempTSTSnum],"??") != NULL || strstr(TSTS_CONNECTION_EQb[tempTSTSnum],"DC") != NULL ){
            TSTS_DC_num[tempTSTSnum] = 1;
            if( strstr(TSTS_CONNECTION_EQa[tempTSTSnum],"??") != NULL || strstr(TSTS_CONNECTION_EQa[tempTSTSnum],"DC") != NULL ){
              TSTS_CONNECTION_EQa_num[tempTSTSnum] = -1;
            }else{
              TSTS_CONNECTION_EQa_num[tempTSTSnum] = atoi(TSTS_CONNECTION_EQa[tempTSTSnum]);
            }
            if( strstr(TSTS_CONNECTION_EQb[tempTSTSnum],"??") != NULL || strstr(TSTS_CONNECTION_EQb[tempTSTSnum],"DC") != NULL ){
              TSTS_CONNECTION_EQb_num[tempTSTSnum] = -1;
            }else{
              TSTS_CONNECTION_EQb_num[tempTSTSnum] = atoi(TSTS_CONNECTION_EQb[tempTSTSnum]);
            }
          }else{
            TSTS_CONNECTION_EQa_num[tempTSTSnum] = atoi(TSTS_CONNECTION_EQa[tempTSTSnum]);
            TSTS_CONNECTION_EQb_num[tempTSTSnum] = atoi(TSTS_CONNECTION_EQb[tempTSTSnum]);
          }
          TSTS_composition_EQa_num[tempTSTSnum] = true_compositionnum[TSTS_CONNECTION_EQa_num[tempTSTSnum]];
          TSTS_composition_EQb_num[tempTSTSnum] = true_compositionnum[TSTS_CONNECTION_EQb_num[tempTSTSnum]];
          if( TSTS_composition_EQa_num[tempTSTSnum] == TSTS_composition_EQb_num[tempTSTSnum] && TSTS_DC_num[tempTSTSnum] == 0){
            TSTS_same_composition_list[tempTSTSnum] = 1;
          }
          tempTSTSnum++;
        }
      }
      for(j = 1;j < totalcompositionnum+1;j++){ // "j" runs "true_compositionnum"
        fprintf(output,"list of TS connected from composition type %d ",j);
        if( no_babel_mode == 0 ){
          fprintf(output,"( Iseq ");
          tempnum = -2;
          i = 0;
            for(k = 0;k < totalEQnum;k++){
              if( true_compositionnum[k] == j && tempnum < inchi_composition_groupnum[k] ){
                fprintf(output,"%d ",inchi_composition_groupnum[k]+1);
                temp_babel_sequence_num_list[i] = inchi_composition_groupnum[k];
                i++;
                tempnum = inchi_composition_groupnum[k];
              }
            }
          temp_babel_sequence_num_listnum = i;
          i = 0;
          fprintf(output,")\n");
        }else{
          fprintf(output,"\n");
        }
        for(k = 0;k < totalTSTSnum;k++){
          if( TSTS_composition_EQa_num[k] == j && TSTS_same_composition_list[k] == 0 && TSTS_DC_num[k] == 0 ){
            if( no_babel_mode == 0 ){
              fprintf(output," TS %*d\t: %*d - %*d\t( -> type %*d : Iseq %*d )",totalTSTS_digitnum,k,totalEQ_digitnum,TSTS_CONNECTION_EQa_num[k],totalEQ_digitnum,TSTS_CONNECTION_EQb_num[k],totalcomposition_digitnum,TSTS_composition_EQb_num[k],totalcomposition_digitnum,inchi_composition_groupnum[TSTS_CONNECTION_EQb_num[k]]+1);
              if( TSlist2_mode == 0 ){ // added 9507
                    fprintf(output," : E = %17.12lf (dE = %.1lf kJ/mol)",TSTS_Bare_Energy[k], 2625.49962 * (TSTS_Bare_Energy[k] - Bare_Energy[TSTS_CONNECTION_EQa_num[k]]));
              }else if( TSlist2_mode == 1 ){
                  fprintf(output," : Rel.G = %.1lf (dG = %.1lf kJ/mol)", TS_RT_Gibbs_Energy[k], TS_RT_Gibbs_Energy[k] - RT_Gibbs_Energy[TSTS_CONNECTION_EQa_num[k]]);
              } //####################################
              for(m = 0;m < temp_babel_sequence_num_listnum;m++){
                if( inchi_composition_groupnum[TSTS_CONNECTION_EQb_num[k]] == temp_babel_sequence_num_list[m] ){
                  fprintf(output," : the same Iseq connection");
                  break;
                }
              }
              fprintf(output,"\n");
            }else{
              fprintf(output," TS %*d\t: %*d - %*d\t( -> type %*d )",totalTSTS_digitnum,k,totalEQ_digitnum,TSTS_CONNECTION_EQa_num[k],totalEQ_digitnum,TSTS_CONNECTION_EQb_num[k],totalcomposition_digitnum,TSTS_composition_EQb_num[k]);
              if( TSlist2_mode == 0 ){ // added 9507
                fprintf(output," : E = %17.12lf (dE = %.1lf kJ/mol)\n",TSTS_Bare_Energy[k], 2625.49962 * (TSTS_Bare_Energy[k] - Bare_Energy[TSTS_CONNECTION_EQa_num[k]]));
              }else if( TSlist2_mode == 1 ){
                fprintf(output," : Rel.G = %.1lf (dG = %.1lf kJ/mol)\n", TS_RT_Gibbs_Energy[k], TS_RT_Gibbs_Energy[k] - RT_Gibbs_Energy[TSTS_CONNECTION_EQa_num[k]]);
              }
            }
          }else if( TSTS_composition_EQb_num[k] == j && TSTS_same_composition_list[k] == 0 && TSTS_DC_num[k] == 0 ){
            if( no_babel_mode == 0 ){
              fprintf(output," TS %*d\t: %*d - %*d\t( -> type %*d : Iseq %*d )",totalTSTS_digitnum,k,totalEQ_digitnum,TSTS_CONNECTION_EQb_num[k],totalEQ_digitnum,TSTS_CONNECTION_EQa_num[k],totalcomposition_digitnum,TSTS_composition_EQa_num[k],totalcomposition_digitnum,inchi_composition_groupnum[TSTS_CONNECTION_EQa_num[k]]+1);
              if( TSlist2_mode == 0 ){ // added 9507
                  fprintf(output," : E = %17.12lf (dE = %.1lf kJ/mol)",TSTS_Bare_Energy[k], 2625.49962 * (TSTS_Bare_Energy[k] - Bare_Energy[TSTS_CONNECTION_EQb_num[k]]));
              }else if( TSlist2_mode == 1 ){
                  fprintf(output," : Rel.G = %.1lf (dG = %.1lf kJ/mol)", TS_RT_Gibbs_Energy[k], TS_RT_Gibbs_Energy[k] - RT_Gibbs_Energy[TSTS_CONNECTION_EQb_num[k]]);
              } //####################################
              for(m = 0;m < temp_babel_sequence_num_listnum;m++){
                if( inchi_composition_groupnum[TSTS_CONNECTION_EQa_num[k]] == temp_babel_sequence_num_list[m] ){
                  fprintf(output," : the same Iseq connection");
                  break;
                }
              }
              fprintf(output,"\n");
            }else{
              fprintf(output," TS %*d\t: %*d - %*d\t( -> type %*d )",totalTSTS_digitnum,k,totalEQ_digitnum,TSTS_CONNECTION_EQb_num[k],totalEQ_digitnum,TSTS_CONNECTION_EQa_num[k],totalcomposition_digitnum,TSTS_composition_EQa_num[k]);
              if( TSlist2_mode == 0 ){ // added 9507
                fprintf(output," : E = %17.12lf (dE = %.1lf kJ/mol)\n",TSTS_Bare_Energy[k], 2625.49962 * (TSTS_Bare_Energy[k] - Bare_Energy[TSTS_CONNECTION_EQb_num[k]]));
              }else if( TSlist2_mode == 1 ){
                fprintf(output," : Rel.G = %.1lf (dG = %.1lf kJ/mol)\n", TS_RT_Gibbs_Energy[k], TS_RT_Gibbs_Energy[k] - RT_Gibbs_Energy[TSTS_CONNECTION_EQb_num[k]]);
              }
            }
          }
          i = 0;
        }
        fprintf(output,"-------------------------------------------------------------------------\n");
      }
      fprintf(output,"\n=========================================================================\nlist of TS connected to the same composition type : (dE = left -> right, right -> left)\n");
      for(j = 0;j < totalTSTSnum;j++){
        if( TSTS_same_composition_list[j] == 1 && TSTS_DC_num[j] == 0 ){
          fprintf(output," TS %*d\t: %*d - %*d\t( in type %*d ) : E = %17.12lf (dE = %.1lf, %.1lf (%.1lf) kJ/mol)\n", totalTSTS_digitnum, j, totalEQ_digitnum, TSTS_CONNECTION_EQa_num[j], totalEQ_digitnum, TSTS_CONNECTION_EQb_num[j], totalcomposition_digitnum, TSTS_composition_EQa_num[j], TSTS_Bare_Energy[j], 2625.49962 * (TSTS_Bare_Energy[j] - Bare_Energy[TSTS_CONNECTION_EQa_num[j]]), 2625.49962 * (TSTS_Bare_Energy[j] - Bare_Energy[TSTS_CONNECTION_EQb_num[j]]), 2625.49962 * (TSTS_Bare_Energy[j] - sorted_Bare_Energy[0]));
        }
      }
      fprintf(output,"=========================================================================\nlist of TS connected to \"??\" or \"DC\" : (dE = left -> right, right -> left)\n");
      for(j = 0;j < totalTSTSnum;j++){
        if( TSTS_DC_num[j] == 1 ){
          // fprintf(output," TS %*d\t: %s - %s\n",totalTSTS_digitnum,j,TSTS_CONNECTION_EQa[j],TSTS_CONNECTION_EQb[j]);
          fprintf(output," TS %*d\t: %s - %s\t: E = %17.12lf (dE = %.1lf, %.1lf (%.1lf) kJ/mol)\n", totalTSTS_digitnum, j, TSTS_CONNECTION_EQa[j],TSTS_CONNECTION_EQb[j], TSTS_Bare_Energy[j], 2625.49962 * (TSTS_Bare_Energy[j] - Bare_Energy[TSTS_CONNECTION_EQa_num[j]]), 2625.49962 * (TSTS_Bare_Energy[j] - Bare_Energy[TSTS_CONNECTION_EQb_num[j]]), 2625.49962 * (TSTS_Bare_Energy[j] - sorted_Bare_Energy[0]));
        }
      }
      fprintf(output,"=========================================================================\n");
      if( TSlist2_mode == 1 ){ // added 9507
        fprintf(output,"\n%s  The values of energy are Gibbs free energy at room temperature.\n  (See \"Relat.Energy ( 298.15 K)\" in %s_TS_list_2.log)\n%s", asterisk, JOBNAME, asterisk);
      }
      fclose(output);
      fclose(fp);





    
      //////////// 9615 added (TS list for Cytoscape) (begin)
      if( cytoscapelistgen_option == 1 ){
        sprintf(cytofilename, "EQCS_%s_cytoTSLOG.csv", JOBNAME);
        cytoscapeTSlist = fopen(cytofilename,"w"); //open the output file
        if(cytoscapeTSlist == NULL) { //failed to open, return NULL
          printf( "failed to open %s!\n", cytofilename);
          return -1;
        }
        fprintf(cytoscapeTSlist, "source,target,energy(ratio)\n");
        for(j = 0;j < totalTSTSnum ;j++){
          fprintf(cytoscapeTSlist,"%s,%s,%10.6lf\n",TSTS_CONNECTION_EQa[j],TSTS_CONNECTION_EQb[j],(TSTS_Bare_Energy[j] - sorted_Bare_Energy[0])/(TSTS_MAXenergy - sorted_Bare_Energy[0]));
        }
        fclose(cytoscapeTSlist);
      }
      //////////// 9615 added (TS list for Cytoscape) (end)
      
    }
  } // if TS_list exists <<*end*>>
  
  ////////////////////////////////////////  vvv PT vvv  ////////////////////////////////////////

  sprintf(PT_fname,"%s_PT_list.log",JOBNAME);
  fp = fopen(PT_fname,"r"); // open TS_list to read
  if(fp != NULL) { // if PT_list exists <<*begin*>>
    PTTS_on = 1;
    while(fgets(str, N, fp) != NULL){
      if( strstr( str , "CONNECTION" ) != NULL ){
        tempPTTSnum++;
      }
    }
    totalPTTSnum = tempPTTSnum;
    tempPTTSnum = 0; /// initialize "tempPTTSnum"
    fclose(fp);

    totalPTTS_digitnum = digit(totalPTTSnum);

    if( totalPTTSnum > 0 ){
      printf("\r  runing... (PT-sorting)");
      char PTTS_CONNECTION_EQa[totalPTTSnum][N],PTTS_CONNECTION_EQb[totalPTTSnum][N];
      double PTTS_Forced_Energy[totalPTTSnum], PTTS_Bare_Energy[totalPTTSnum], PTTS_Zero_Energy[totalPTTSnum];
      int PTTS_CONNECTION_EQa_num[totalPTTSnum],PTTS_CONNECTION_EQb_num[totalPTTSnum];
      int PTTS_composition_EQa_num[totalPTTSnum],PTTS_composition_EQb_num[totalPTTSnum];
      int PTTS_DC_num[totalPTTSnum]; // 0 (not to DC) or 1 (connected to DC)
      int PTTS_same_composition_list[totalPTTSnum]; // 0 (not same) or 1 (connected to the same composition)
      for(j = 0;j < totalPTTSnum;j++){
        PTTS_same_composition_list[j] = 0;
        PTTS_DC_num[j] = 0;
      }
      sprintf(tempPT_output_fname, "tempEQcomposition_PTlist_%s.log", JOBNAME);
      output = fopen(tempPT_output_fname,"w"); // open the output file
      if(output == NULL) { // failed to open, return NULL
        printf("failed to open %s!\n", tempPT_output_fname);
        return -1;
      }
      if( no_babel_mode == 0 ){
        fprintf(output,"*************************************************************************\n   EQ composition PT list (version %s)\n      \"%s\"\n*************************************************************************\n                                            produced by H.N. (%s)\n   the number of EQ            : %d\n   the number of atom          : %d\n   the number of element       : %d\n   the number of composition   : %d\n   the number of sequence type : %d\n\n-------------------------------------------------------------------------\n\n   << list of compositions >>\n\n", software_version, JOBNAME, developed_date, totalEQnum, atoms, totalelementnamenum, totalcompositionnum, total_inchi_compositionnum);
      }else{
        fprintf(output,"*************************************************************************\n   EQ composition PT list (version %s)\n      \"%s\"\n*************************************************************************\n                                            produced by H.N. (%s)\n   the number of EQ            : %d\n   the number of atom          : %d\n   the number of element       : %d\n   the number of composition   : %d\n\n\n-------------------------------------------------------------------------\n\n   << list of compositions >>\n\n", software_version, JOBNAME, developed_date, totalEQnum, atoms, totalelementnamenum, totalcompositionnum);
      }

      for(j = 0;j < totalcompositionnum;j++){ // list up numbers of EQ whose composition type is "j"
        fprintf(output,"   type %d EQ\t:",j+1);
        i = 0;
        for(k = 0;k < totalEQnum;k++){
          if( compositionnum[k] == j ){
            fprintf(output," %d",k);
            i++;
          }
        }
        fprintf(output," (%d EQ(s))\n             \t:",i);
        i = 0;
        tempnum=0;
        for(k = 0;k < bonding_groups_num[j];k++){
          if( tempcomposition_list[j][k] == NULL ){
            break;
          }else{
            fprintf(output," %s", tempcomposition_list[j][k]);
            tempnum++;
          }
        }
        fprintf(output," (%d molecule(s))\n",tempnum);
      }
      tempnum=0;
      fprintf(output,"\n"); // added 8c31
      if( no_babel_mode == 0 ){
        fprintf(output,"-------------------------------------------------------------------------\n\n   << list of sequences >>\n\n");
        for(j = 0;j < total_inchi_compositionnum;j++){ // list up numbers of EQ whose sequence type is "j" // added 8c31
          fprintf(output,"   Iseq %*d EQ\t:",totalcomposition_digitnum,j+1);
          for(k = 0;k < totalEQnum;k++){
            if( inchi_composition_groupnum[k] == j ){
              fprintf(output," %d",k);
              i++;
            }
          }
          fprintf(output," (%d EQ(s))\n",i);
          i = 0;
        }
      }
      fprintf(output,"\n=========================================================================\n\n");

      fp = fopen(PT_fname,"r");
      while(fgets(str, N, fp) != NULL){
        if( strstr( str , "Energy" ) != NULL ){
          if( strstr( str , ":" ) != NULL ){
            sscanf(str, "Energy    = %lf (%lf : %lf)\n", &PTTS_Forced_Energy[tempPTTSnum], &PTTS_Bare_Energy[tempPTTSnum], &PTTS_Zero_Energy[tempPTTSnum]); // scan TS barrier
            if( PTTS_Bare_Energy[tempPTTSnum] > -0.1 ){
              PTTS_Bare_Energy[tempPTTSnum] = PTTS_Forced_Energy[tempPTTSnum];
            }
          }else{
            sscanf(str, "Energy    = %lf\n", &PTTS_Bare_Energy[tempPTTSnum]); // scan TS barrier
          }
          if( PTTS_MAXenergy < PTTS_Bare_Energy[tempPTTSnum] ){
            PTTS_MAXenergy = PTTS_Bare_Energy[tempPTTSnum];
          }
        }
        if( strstr( str , "CONNECTION" ) != NULL ){
          sscanf(str, "%s %s %s %s %s\n",temp_char[0],temp_char[1],PTTS_CONNECTION_EQa[tempPTTSnum],temp_char[2],PTTS_CONNECTION_EQb[tempPTTSnum]); // scan "CONNECTION : 0 - ??"
          if( strstr(PTTS_CONNECTION_EQa[tempPTTSnum],"??") != NULL || strstr(PTTS_CONNECTION_EQa[tempPTTSnum],"DC") != NULL || strstr(PTTS_CONNECTION_EQb[tempPTTSnum],"??") != NULL || strstr(PTTS_CONNECTION_EQb[tempPTTSnum],"DC") != NULL ){
            PTTS_DC_num[tempPTTSnum] = 1;
            if( strstr(PTTS_CONNECTION_EQa[tempPTTSnum],"??") != NULL || strstr(PTTS_CONNECTION_EQa[tempPTTSnum],"DC") != NULL ){
              PTTS_CONNECTION_EQa_num[tempPTTSnum] = -1;
            }else{
              PTTS_CONNECTION_EQa_num[tempPTTSnum] = atoi(PTTS_CONNECTION_EQa[tempPTTSnum]);
            }
            if( strstr(PTTS_CONNECTION_EQb[tempPTTSnum],"??") != NULL || strstr(PTTS_CONNECTION_EQa[tempPTTSnum],"DC") != NULL ){
              PTTS_CONNECTION_EQb_num[tempPTTSnum] = -1;
            }else{
              PTTS_CONNECTION_EQb_num[tempPTTSnum] = atoi(PTTS_CONNECTION_EQb[tempPTTSnum]);
            }
          }else{
            PTTS_CONNECTION_EQa_num[tempPTTSnum] = atoi(PTTS_CONNECTION_EQa[tempPTTSnum]);
            PTTS_CONNECTION_EQb_num[tempPTTSnum] = atoi(PTTS_CONNECTION_EQb[tempPTTSnum]);
          }
          PTTS_composition_EQa_num[tempPTTSnum] = true_compositionnum[PTTS_CONNECTION_EQa_num[tempPTTSnum]];
          PTTS_composition_EQb_num[tempPTTSnum] = true_compositionnum[PTTS_CONNECTION_EQb_num[tempPTTSnum]];
          if( PTTS_composition_EQa_num[tempPTTSnum] == PTTS_composition_EQb_num[tempPTTSnum] && PTTS_DC_num[tempPTTSnum] == 0){
            PTTS_same_composition_list[tempPTTSnum] = 1;
          }
          tempPTTSnum++;
        }
      }
      for(j = 1;j < totalcompositionnum+1;j++){ // "j" runs "true_compositionnum"
        fprintf(output,"list of PT connected from composition type %d ",j);
        if( no_babel_mode == 0 ){
          fprintf(output,"( Iseq ");
          tempnum = -2;
          i = 0;
            for(k = 0;k < totalEQnum;k++){
              if( true_compositionnum[k] == j && tempnum < inchi_composition_groupnum[k] ){
                fprintf(output,"%d ",inchi_composition_groupnum[k]+1);
                temp_babel_sequence_num_list[i] = inchi_composition_groupnum[k];
                i++;
                tempnum = inchi_composition_groupnum[k];
              }
            }
          temp_babel_sequence_num_listnum = i;
          i = 0;
          fprintf(output,")\n");
        }else{
          fprintf(output,"\n");
        }
        for(k = 0;k < totalPTTSnum;k++){ // "k" runs as a tempPTTSnum
          if( PTTS_composition_EQa_num[k] == j && PTTS_same_composition_list[k] == 0 && PTTS_DC_num[k] == 0 ){
            if( no_babel_mode == 0 ){
              fprintf(output," PT %*d\t: %*d - %*d\t( -> type %*d : Iseq %*d )",totalPTTS_digitnum,k,totalEQ_digitnum,PTTS_CONNECTION_EQa_num[k],totalEQ_digitnum,PTTS_CONNECTION_EQb_num[k],totalcomposition_digitnum,PTTS_composition_EQb_num[k],totalcomposition_digitnum,inchi_composition_groupnum[PTTS_CONNECTION_EQb_num[k]]+1);
              fprintf(output," : E = %17.12lf (dE = %3.1lf, dErev = %3.1lf (dErel = %3.1lf) kJ/mol)",PTTS_Bare_Energy[k], 2625.49962 * (PTTS_Bare_Energy[k] - Bare_Energy[PTTS_CONNECTION_EQa_num[k]]), 2625.49962 * (PTTS_Bare_Energy[k] - Bare_Energy[PTTS_CONNECTION_EQb_num[k]]), 2625.49962 * (PTTS_Bare_Energy[k] - sorted_Bare_Energy[0]));
              for(m = 0;m < temp_babel_sequence_num_listnum;m++){
                if( inchi_composition_groupnum[PTTS_CONNECTION_EQb_num[k]] == temp_babel_sequence_num_list[m] ){
                  fprintf(output," : the same Iseq connection");
                  break;
                }
              }
              fprintf(output,"\n");
            }else{
              fprintf(output," PT %*d\t: %*d - %*d\t( -> type %*d ) : E = %17.12lf (dE = %3.1lf, dErev = %3.1lf (dErel = %3.1lf) kJ/mol)\n", totalPTTS_digitnum, k, totalEQ_digitnum, PTTS_CONNECTION_EQa_num[k], totalEQ_digitnum, PTTS_CONNECTION_EQb_num[k], totalcomposition_digitnum, PTTS_composition_EQb_num[k], totalcomposition_digitnum, PTTS_composition_EQb_num[k], PTTS_Bare_Energy[k], 2625.49962 * (PTTS_Bare_Energy[k] - Bare_Energy[PTTS_CONNECTION_EQa_num[k]]), 2625.49962 * (PTTS_Bare_Energy[k] - Bare_Energy[PTTS_CONNECTION_EQb_num[k]]), 2625.49962 * (PTTS_Bare_Energy[k] - sorted_Bare_Energy[0]));
            }
          }else if( PTTS_composition_EQb_num[k] == j && PTTS_same_composition_list[k] == 0 && PTTS_DC_num[k] == 0 ){
            if( no_babel_mode == 0 ){
              fprintf(output," PT %*d\t: %*d - %*d\t( -> type %*d : Iseq %*d )",totalPTTS_digitnum,k,totalEQ_digitnum,PTTS_CONNECTION_EQb_num[k],totalEQ_digitnum,PTTS_CONNECTION_EQa_num[k],totalcomposition_digitnum,PTTS_composition_EQa_num[k], totalcomposition_digitnum, inchi_composition_groupnum[PTTS_CONNECTION_EQa_num[k]]+1);
              fprintf(output," : E = %17.12lf (dE = %3.1lf, dErev = %3.1lf (dErel = %3.1lf) kJ/mol)", PTTS_Bare_Energy[k], 2625.49962 * (PTTS_Bare_Energy[k] - Bare_Energy[PTTS_CONNECTION_EQb_num[k]]), 2625.49962 * (PTTS_Bare_Energy[k] - Bare_Energy[PTTS_CONNECTION_EQa_num[k]]), 2625.49962 * (PTTS_Bare_Energy[k] - sorted_Bare_Energy[0]));
              for(m = 0;m < temp_babel_sequence_num_listnum;m++){
                if( inchi_composition_groupnum[PTTS_CONNECTION_EQa_num[k]] == temp_babel_sequence_num_list[m] ){
                  fprintf(output," : the same Iseq connection");
                  break;
                }
              }
              fprintf(output,"\n");
            }else{
              fprintf(output," PT %*d\t: %*d - %*d\t( -> type %*d ) : E = %17.12lf (dE = %3.1lf, dErev = %3.1lf (dErel = %3.1lf) kJ/mol)\n", totalPTTS_digitnum, k, totalEQ_digitnum, PTTS_CONNECTION_EQb_num[k], totalEQ_digitnum, PTTS_CONNECTION_EQa_num[k], totalcomposition_digitnum, PTTS_composition_EQa_num[k], PTTS_Bare_Energy[k], 2625.49962 * (PTTS_Bare_Energy[k] - Bare_Energy[PTTS_CONNECTION_EQb_num[k]]), 2625.49962 * (PTTS_Bare_Energy[k] - Bare_Energy[PTTS_CONNECTION_EQa_num[k]]), 2625.49962 * (PTTS_Bare_Energy[k] - sorted_Bare_Energy[0]));
            }
          }
          i = 0;
        }
        fprintf(output,"-------------------------------------------------------------------------\n");
      }
      fprintf(output,"\n=========================================================================\nlist of PT connected to the same composition type : (dE = left -> right, right -> left)\n");
      for(j = 0;j < totalPTTSnum;j++){
        if( PTTS_same_composition_list[j] == 1 && PTTS_DC_num[j] == 0 ){
          fprintf(output," PT %*d\t: %*d - %*d\t( in type %*d ) : E = %17.12lf (dE = %3.1lf, %3.1lf (dErel = %3.1lf) kJ/mol)\n", totalPTTS_digitnum, j, totalEQ_digitnum, PTTS_CONNECTION_EQa_num[j], totalEQ_digitnum, PTTS_CONNECTION_EQb_num[j], totalcomposition_digitnum, PTTS_composition_EQa_num[j], PTTS_Bare_Energy[j], 2625.49962 * (PTTS_Bare_Energy[j] - Bare_Energy[PTTS_CONNECTION_EQa_num[j]]), 2625.49962 * (PTTS_Bare_Energy[j] - Bare_Energy[PTTS_CONNECTION_EQb_num[j]]), 2625.49962 * (PTTS_Bare_Energy[j] - sorted_Bare_Energy[0]));
        }
      }
      fprintf(output,"=========================================================================\nlist of PT connected to \"??\" or \"DC\" : (dE = left -> right, right -> left)\n");
      for(j = 0;j < totalPTTSnum;j++){
        if( PTTS_DC_num[j] == 1 ){
          // fprintf(output," PT %*d\t: %s - %s\n",totalPTTS_digitnum,j,PTTS_CONNECTION_EQa[j],PTTS_CONNECTION_EQb[j]);
          fprintf(output," PT %*d\t: %s - %s\t: E = %17.12lf (dE = %3.1lf, %3.1lf (dErel = %3.1lf) kJ/mol)\n", totalPTTS_digitnum, j, PTTS_CONNECTION_EQa[j],PTTS_CONNECTION_EQb[j], PTTS_Bare_Energy[j], 2625.49962 * (PTTS_Bare_Energy[j] - Bare_Energy[PTTS_CONNECTION_EQa_num[j]]), 2625.49962 * (PTTS_Bare_Energy[j] - Bare_Energy[PTTS_CONNECTION_EQb_num[j]]), 2625.49962 * (PTTS_Bare_Energy[j] - sorted_Bare_Energy[0]));
        }
      }
      fprintf(output,"=========================================================================\n");
      fclose(output);
      fclose(fp);
      
      //////////// 9615 added (PT list for Cytoscape) (begin)
      if( cytoscapelistgen_option == 1 ){
        sprintf(cytofilename, "EQCS_%s_cytoPTLOG.csv", JOBNAME);
        cytoscapePTlist = fopen(cytofilename,"w"); //open the output file
        if(cytoscapePTlist == NULL) { //failed to open, return NULL
          printf( "failed to open %s!\n", cytofilename);
          return -1;
        }
        fprintf(cytoscapePTlist, "source,target,energy(ratio)\n");
        for(j = 0;j < totalPTTSnum ;j++){
          fprintf(cytoscapePTlist,"%s,%s,%10.6lf\n",PTTS_CONNECTION_EQa[j],PTTS_CONNECTION_EQb[j],(PTTS_Bare_Energy[j] - sorted_Bare_Energy[0])/(PTTS_MAXenergy - sorted_Bare_Energy[0]));
        }
        fclose(cytoscapePTlist);
      }
      //////////// 9615 added (PT list for Cytoscape) (end)
      
    }
  } // if PT_list exists <<*end*>>
  
    sprintf(tempTS_output_fname, "tempEQcomposition_TSlist_%s.log", JOBNAME);
    fp = fopen(tempTS_output_fname,"r"); // open the output file
    if(fp != NULL) { // failed to open, return NULL
      line_num = 1;
      if( TSlist2_mode == 1 ){ // added 9507
        sprintf(TS_output_fname, "EQCS_%s_TSLOG2.log", JOBNAME);
      }else{
        sprintf(TS_output_fname, "EQCS_%s_TSLOG.log", JOBNAME);
      }
      re_output = fopen(TS_output_fname,"w"); // open the output file
      while(fgets(str, N, fp) != NULL){
        fprintf(re_output,"%s",str);
        if( line_num == 10 ){
          fprintf(re_output,"   the number of TS            : %d",totalTSTSnum);
          if( TSTS_on == 0 ){
            fprintf(re_output,"\t(TS_list.log doesn't exist)");
          }
          fprintf(re_output,"\n   the number of PT            : %d",totalPTTSnum);
          if( PTTS_on == 0 ){
            fprintf(re_output,"\t(PT_list.log doesn't exist)");
          }
          fprintf(re_output,"\n");
        }
        line_num++;
      }
      fclose(re_output);
      fclose(fp);
      remove(tempTS_output_fname);
    }

    sprintf(tempPT_output_fname, "tempEQcomposition_PTlist_%s.log", JOBNAME);
    fp = fopen(tempPT_output_fname,"r"); // open the output file
    if(fp != NULL) { // failed to open, return NULL
      line_num = 1;
      sprintf(PT_output_fname, "EQCS_%s_PTLOG.log", JOBNAME);
      re_output = fopen(PT_output_fname,"w"); // open the output file
      while(fgets(str, N, fp) != NULL){
        fprintf(re_output,"%s",str);
        if( line_num == 10 ){
          fprintf(re_output,"   the number of TS            : %d",totalTSTSnum);
          if( TSTS_on == 0 ){
            fprintf(re_output,"\t(TS_list.log doesn't exist)");
          }
          fprintf(re_output,"\n   the number of PT            : %d",totalPTTSnum);
          if( PTTS_on == 0 ){
            fprintf(re_output,"\t(PT_list.log doesn't exist)");
          }
          fprintf(re_output,"\n");
        }
        line_num++;
      }
      fclose(re_output);
      fclose(fp);
      remove(tempPT_output_fname);
    }
  
  // Standard output
  if( TSTS_on != 0 ){
    printf("\r   the number of TS : %d / TS_MAXenergy = %17.12lf\n",totalTSTSnum,TSTS_MAXenergy);
  }
  if( PTTS_on != 0 ){
    printf("\r   the number of PT : %d / PT_MAXenergy = %17.12lf\n",totalPTTSnum,PTTS_MAXenergy);
  }
  fprintf(re_output,"\n");
  if( TSTS_on != 0 && PTTS_on != 0 ){
    printf("   TS/PT (Energy normalization ratio) = %7.6lf\n",(TSTS_MAXenergy - sorted_Bare_Energy[0])/(PTTS_MAXenergy - sorted_Bare_Energy[0]));
  }
  printf("\r  *EQ composition sorting was completed. ");
  

    ///// energy_sort_mode part (begin)
    if(energy_sort_mode == 1 && EQlist2_mode == 0 ){
      if( density_sort_mode == 1 && no_babel_mode == 0 ){
        fprintf(energy_sort,"rank \tEQ# \tISeq# \trel. E(el) (kJ/mol) ( hartree : g/cm^3 )\n");
      }else if( no_babel_mode == 0 ){
        fprintf(energy_sort,"rank \tEQ# \tISeq# \trel. E(el) (kJ/mol) ( hartree )\n");
      }else{
        fprintf(energy_sort,"rank \tEQ# \ttype# \trel. E(el) (kJ/mol) ( hartree )\n");
      }
      for(k = 0;k < totalEQnum ;k++){
        for(m = 0;m < totalEQnum ;m++){
          if( Energy_rank_EQnum[m] == k ){
            if( density_sort_mode == 1 && no_babel_mode == 0 ){
              fprintf(energy_sort," %*d \tEQ %-*d \tISeq %-*d \t%12.7lf ( %17.12lf : %6.3lf )\n", totalEQ_digitnum, k+1, totalEQ_digitnum, Energy_rank_EQnum[k], totalcomposition_digitnum, inchi_composition_groupnum[Energy_rank_EQnum[k]]+1, (sorted_Bare_Energy[k] - sorted_Bare_Energy[0]) * 2625.49962, sorted_Bare_Energy[k], cell_density[Energy_rank_EQnum[k]]); // (9419 added / 9613 modified)
            }else if( no_babel_mode == 0 ){
              fprintf(energy_sort," %*d \tEQ %-*d \tISeq %-*d \t%12.7lf ( %17.12lf )\n", totalEQ_digitnum, k+1, totalEQ_digitnum, Energy_rank_EQnum[k], totalcomposition_digitnum, inchi_composition_groupnum[Energy_rank_EQnum[k]]+1, (sorted_Bare_Energy[k] - sorted_Bare_Energy[0]) * 2625.49962, sorted_Bare_Energy[k]); // (9613 modified)
            }else{
              fprintf(energy_sort," %*d \tEQ %-*d \ttype %-*d \t%12.7lf ( %17.12lf )\n", totalEQ_digitnum, k+1, totalEQ_digitnum, Energy_rank_EQnum[k], totalcomposition_digitnum, inchi_composition_groupnum[Energy_rank_EQnum[k]]+1, (sorted_Bare_Energy[k] - sorted_Bare_Energy[0]) * 2625.49962, sorted_Bare_Energy[k]); // (9613 modified)
            }
            // 23606 added
            if(crystal_surface_mode == 0){
              fprintf(energy_sort_xyz,"%d\n# EQ %d ; %12.7lf ( %17.12lf )\n", atoms, Energy_rank_EQnum[k], (sorted_Bare_Energy[k] - sorted_Bare_Energy[0]) * 2625.49962, sorted_Bare_Energy[k]);
              for(j = 0;j < atoms ;j++){
                fprintf(energy_sort_xyz,"%s    %lf    %lf    %lf\n", atomname[Energy_rank_EQnum[k]][j], temp_coordinate_x[Energy_rank_EQnum[k]][j], temp_coordinate_y[Energy_rank_EQnum[k]][j], temp_coordinate_z[Energy_rank_EQnum[k]][j]);
              }
            }else if(crystal_surface_mode == 1){
              fprintf(energy_sort_xyz,"%d\n# EQ %d ; %12.7lf ( %17.12lf )\n", atoms + atoms * (framescale_x * framescale_y * framescale_z - 1), Energy_rank_EQnum[k], (sorted_Bare_Energy[k] - sorted_Bare_Energy[0]) * 2625.49962, sorted_Bare_Energy[k]);
              for(jj = 0;jj < meta_atoms;jj++){
                if( jj != TVatom_linenum[0] && jj != TVatom_linenum[1] && jj != TVatom_linenum[2]){
                  double tempmeta_X = temp_coordinate_x[Energy_rank_EQnum[k]][jj];
                  double tempmeta_Y = temp_coordinate_y[Energy_rank_EQnum[k]][jj];
                  double tempmeta_Z = temp_coordinate_z[Energy_rank_EQnum[k]][jj];
                  fprintf(energy_sort_xyz, "%-2s\t%17.12lf %17.12lf %17.12lf\n",atomname[Energy_rank_EQnum[k]][jj], tempmeta_X, tempmeta_Y, tempmeta_Z);
                }
              }
              for(kk = 0;kk <= framescale_x;kk++){
                for(ll = 0;ll <= framescale_y;ll++){
                  for(mm = 0;mm <= framescale_z;mm++){
                    for(jj = 0;jj < meta_atoms;jj++){
                      if( jj != TVatom_linenum[0] && jj != TVatom_linenum[1] && jj != TVatom_linenum[2] && kk+ll+mm != 0 ){
                        double tempmeta_X = temp_coordinate_x[Energy_rank_EQnum[k]][jj] + kk*temp_coordinate_x[Energy_rank_EQnum[k]][TVatom_linenum[0]] + ll*temp_coordinate_x[Energy_rank_EQnum[k]][TVatom_linenum[1]] + mm*temp_coordinate_x[Energy_rank_EQnum[k]][TVatom_linenum[2]];
                        double tempmeta_Y = temp_coordinate_y[Energy_rank_EQnum[k]][jj] + kk*temp_coordinate_y[Energy_rank_EQnum[k]][TVatom_linenum[0]] + ll*temp_coordinate_y[Energy_rank_EQnum[k]][TVatom_linenum[1]] + mm*temp_coordinate_y[Energy_rank_EQnum[k]][TVatom_linenum[2]];
                        double tempmeta_Z = temp_coordinate_z[Energy_rank_EQnum[k]][jj] + kk*temp_coordinate_z[Energy_rank_EQnum[k]][TVatom_linenum[0]] + ll*temp_coordinate_z[Energy_rank_EQnum[k]][TVatom_linenum[1]] + mm*temp_coordinate_z[Energy_rank_EQnum[k]][TVatom_linenum[2]];
                        fprintf(energy_sort_xyz, "%-2s\t%17.12lf %17.12lf %17.12lf\n",atomname[Energy_rank_EQnum[k]][jj], tempmeta_X, tempmeta_Y, tempmeta_Z);
                      }
                    }
                  }
                }
              }
            }
            fprintf(energy_sort_xyz,"\n");
            // printf("\n%d", k);
            break;
          }
        }
      }
      fclose(energy_sort);
      fclose(energy_sort_xyz);
    }
    ///// energy_sort_mode part (end)
    
  } ///// ##### composition sorting mode is on (end)

  printf("(EQCS ver %s)\n", software_version);

  return 0;
}


////////// << function 1 "strcount" >> ////////// search *str2 in *str1
int strcount(char *str1, char *str2) {
  char *start_ptr;  /* search begin point */
  char *search_ptr; /* search end point */
  int occurrences_number = 0;

  start_ptr = str1;
  search_ptr = strstr(start_ptr, str2);
  while( search_ptr != NULL ){
    start_ptr = search_ptr + strlen(str2);
    occurrences_number++;
    search_ptr = strstr(start_ptr, str2);
  }
  return occurrences_number;
}

////////// << function 2 "mystrcmp" >> ////////// compare *str2 and *str1
int mystrcmp(char *str1, char *str2){
  int i;
  for(i = 0;*(str1 + i) == *(str2 + i);i++){
    if(*(str1 + i) == '\0'){
      return 0;
    }
  }
  return *(str1 + i) - *(str2 + i);
} 

////////// << function 3 "digit" >> ////////// determinate a number of digit
int digit(int number){
  int di = 0;
  int tempval = number;
  while( tempval != 0 ){
    tempval = tempval / 10;
    di++;
  }
  return di;
}

////////// << function 4 "lntrim" >> ////////// display progress ( number should be a integer from 0 to 100 )
void lntrim(char *str){
  int i = 0;
  while(1){
    if(str[i] == '\n'){
      str[i] = '\0';
      break;
    }
    i++;
  }
}

////////// << function 5 "radius" >> ////////// determinate the two atoms have a Covalent bond
double radius(char *a){
  if(strcmp(a,"HX") == 0){
//    return 0.31; // original value
    return 0.32;
  }
  if(strcmp(a,"HeX") == 0){
    return 0.28;
  }
  if(strcmp(a,"LiX") == 0){
    return 1.28;
  }
  if(strcmp(a,"BeX") == 0){
    return 0.96;
  }
  if(strcmp(a,"BX") == 0){
    return 0.84;
  }
  if(strcmp(a,"CX") == 0){
    return 0.76;
  }
  if(strcmp(a,"NX") == 0){
    return 0.71;
  }
  if(strcmp(a,"OX") == 0){
    return 0.66;
  }
  if(strcmp(a,"FX") == 0){
    return 0.57;
  }
  if(strcmp(a,"NeX") == 0){
    return 0.58;
  }
  if(strcmp(a,"NaX") == 0){
    return 1.66;
  }
  if(strcmp(a,"MgX") == 0){
    return 1.41;
  }
  if(strcmp(a,"AlX") == 0){
    return 1.21;
  }
  if(strcmp(a,"SiX") == 0){
//    return 1.11; // original value
    return 1.16;
  }
  if(strcmp(a,"PX") == 0){
    return 1.07;
  }
  if(strcmp(a,"SX") == 0){
    return 1.05;
  }
  if(strcmp(a,"ClX") == 0){
    return 1.02;
  }
  if(strcmp(a,"ArX") == 0){
    return 1.06;
  }
  if(strcmp(a,"KX") == 0){
    return 2.03;
  }
  if(strcmp(a,"CaX") == 0){
    return 1.76;
  }
  if(strcmp(a,"ScX") == 0){
    return 1.7;
  }
  if(strcmp(a,"TiX") == 0){
    return 1.6;
  }
  if(strcmp(a,"VX") == 0){
    return 1.53;
  }
  if(strcmp(a,"CrX") == 0){
    return 1.39;
  }
  if(strcmp(a,"MnX") == 0){
    return 1.39;
  }
  if(strcmp(a,"FeX") == 0){
    return 1.32;
  }
  if(strcmp(a,"CoX") == 0){
    return 1.26;
  }
  if(strcmp(a,"NiX") == 0){
    return 1.24;
  }
  if(strcmp(a,"CuX") == 0){
    return 1.32;
  }
  if(strcmp(a,"ZnX") == 0){
    return 1.22;
  }
  if(strcmp(a,"GaX") == 0){
    return 1.22;
  }
  if(strcmp(a,"GeX") == 0){
    return 1.2;
  }
  if(strcmp(a,"AsX") == 0){
    return 1.19;
  }
  if(strcmp(a,"SeX") == 0){
    return 1.2;
  }
  if(strcmp(a,"BrX") == 0){
    return 1.2;
  }
  if(strcmp(a,"KrX") == 0){
    return 1.16;
  }
  if(strcmp(a,"RbX") == 0){
    return 2.2;
  }
  if(strcmp(a,"SrX") == 0){
    return 1.95;
  }
  if(strcmp(a,"YX") == 0){
    return 1.9;
  }
  if(strcmp(a,"ZrX") == 0){
    return 1.75;
  }
  if(strcmp(a,"NbX") == 0){
    return 1.64;
  }
  if(strcmp(a,"MoX") == 0){
    return 1.54;
  }
  if(strcmp(a,"TcX") == 0){
    return 1.47;
  }
  if(strcmp(a,"RuX") == 0){
    return 1.46;
  }
  if(strcmp(a,"RhX") == 0){
    return 1.42;
  }
  if(strcmp(a,"PdX") == 0){
    return 1.39;
  }
  if(strcmp(a,"AgX") == 0){
    return 1.45;
  }
  if(strcmp(a,"CdX") == 0){
    return 1.44;
  }
  if(strcmp(a,"InX") == 0){
    return 1.42;
  }
  if(strcmp(a,"SnX") == 0){
    return 1.39;
  }
  if(strcmp(a,"SbX") == 0){
    return 1.39;
  }
  if(strcmp(a,"TeX") == 0){
    return 1.38;
  }
  if(strcmp(a,"IX") == 0){
    return 1.39;
  }
  if(strcmp(a,"XeX") == 0){
    return 1.4;
  }
  if(strcmp(a,"CsX") == 0){
    return 2.44;
  }
  if(strcmp(a,"BaX") == 0){
    return 2.15;
  }
  if(strcmp(a,"LaX") == 0){
    return 2.07;
  }
  if(strcmp(a,"CeX") == 0){
    return 2.04;
  }
  if(strcmp(a,"PrX") == 0){
    return 2.03;
  }
  if(strcmp(a,"NdX") == 0){
    return 2.01;
  }
  if(strcmp(a,"PmX") == 0){
    return 1.99;
  }
  if(strcmp(a,"SmX") == 0){
    return 1.98;
  }
  if(strcmp(a,"EuX") == 0){
    return 1.98;
  }
  if(strcmp(a,"GdX") == 0){
    return 1.96;
  }
  if(strcmp(a,"TbX") == 0){
    return 1.94;
  }
  if(strcmp(a,"DyX") == 0){
    return 1.92;
  }
  if(strcmp(a,"HoX") == 0){
    return 1.92;
  }
  if(strcmp(a,"ErX") == 0){
    return 1.89;
  }
  if(strcmp(a,"TmX") == 0){
    return 1.9;
  }
  if(strcmp(a,"YbX") == 0){
    return 1.87;
  }
  if(strcmp(a,"LuX") == 0){
    return 1.87;
  }
  if(strcmp(a,"HfX") == 0){
    return 1.75;
  }
  if(strcmp(a,"TaX") == 0){
    return 1.7;
  }
  if(strcmp(a,"WX") == 0){
    return 1.62;
  }
  if(strcmp(a,"ReX") == 0){
    return 1.51;
  }
  if(strcmp(a,"OsX") == 0){
    return 1.44;
  }
  if(strcmp(a,"IrX") == 0){
    return 1.41;
  }
  if(strcmp(a,"PtX") == 0){
    return 1.36;
  }
  if(strcmp(a,"AuX") == 0){
    return 1.36;
  }
  if(strcmp(a,"HgX") == 0){
    return 1.32;
  }
  if(strcmp(a,"TlX") == 0){
    return 1.45;
  }
  if(strcmp(a,"PbX") == 0){
    return 1.46;
  }
  if(strcmp(a,"BiX") == 0){
    return 1.48;
  }
  if(strcmp(a,"PoX") == 0){
    return 1.4;
  }
  if(strcmp(a,"AtX") == 0){
    return 1.5;
  }
  if(strcmp(a,"RnX") == 0){
    return 1.5;
  }
  if(strcmp(a,"FrX") == 0){
    return 2.6;
  }
  if(strcmp(a,"RaX") == 0){
    return 2.21;
  }
  if(strcmp(a,"AcX") == 0){
    return 2.15;
  }
  if(strcmp(a,"ThX") == 0){
    return 2.06;
  }
  if(strcmp(a,"PaX") == 0){
    return 2.0;
  }
  if(strcmp(a,"UX") == 0){
    return 1.96;
  }
  if(strcmp(a,"NpX") == 0){
    return 1.9;
  }
  if(strcmp(a,"PuX") == 0){
    return 1.87;
  }
  if(strcmp(a,"AmX") == 0){
    return 1.8;
  }
  if(strcmp(a,"CmX") == 0){
    return 1.69;
  }
  return 2.00;
}

////////// << function 6 "Zvalue" >> ////////// return atomic number

int Zvalue(char *a){
  if(strcmp(a,"HX") == 0){
    return 1;
  }
  if(strcmp(a,"HeX") == 0){
    return 2;
  }
  if(strcmp(a,"LiX") == 0){
    return 3;
  }
  if(strcmp(a,"BeX") == 0){
    return 4;
  }
  if(strcmp(a,"BX") == 0){
    return 5;
  }
  if(strcmp(a,"CX") == 0){
    return 6;
  }
  if(strcmp(a,"NX") == 0){
    return 7;
  }
  if(strcmp(a,"OX") == 0){
    return 8;
  }
  if(strcmp(a,"FX") == 0){
    return 9;
  }
  if(strcmp(a,"NeX") == 0){
    return 10;
  }
  if(strcmp(a,"NaX") == 0){
    return 11;
  }
  if(strcmp(a,"MgX") == 0){
    return 12;
  }
  if(strcmp(a,"AlX") == 0){
    return 13;
  }
  if(strcmp(a,"SiX") == 0){
    return 14;
  }
  if(strcmp(a,"PX") == 0){
    return 15;
  }
  if(strcmp(a,"SX") == 0){
    return 16;
  }
  if(strcmp(a,"ClX") == 0){
    return 17;
  }
  if(strcmp(a,"ArX") == 0){
    return 18;
  }
  if(strcmp(a,"KX") == 0){
    return 19;
  }
  if(strcmp(a,"CaX") == 0){
    return 20;
  }
  if(strcmp(a,"ScX") == 0){
    return 21;
  }
  if(strcmp(a,"TiX") == 0){
    return 22;
  }
  if(strcmp(a,"VX") == 0){
    return 23;
  }
  if(strcmp(a,"CrX") == 0){
    return 24;
  }
  if(strcmp(a,"MnX") == 0){
    return 25;
  }
  if(strcmp(a,"FeX") == 0){
    return 26;
  }
  if(strcmp(a,"CoX") == 0){
    return 27;
  }
  if(strcmp(a,"NiX") == 0){
    return 28;
  }
  if(strcmp(a,"CuX") == 0){
    return 29;
  }
  if(strcmp(a,"ZnX") == 0){
    return 30;
  }
  if(strcmp(a,"GaX") == 0){
    return 31;
  }
  if(strcmp(a,"GeX") == 0){
    return 32;
  }
  if(strcmp(a,"AsX") == 0){
    return 33;
  }
  if(strcmp(a,"SeX") == 0){
    return 34;
  }
  if(strcmp(a,"BrX") == 0){
    return 35;
  }
  if(strcmp(a,"KrX") == 0){
    return 36;
  }
  if(strcmp(a,"RbX") == 0){
    return 37;
  }
  if(strcmp(a,"SrX") == 0){
    return 38;
  }
  if(strcmp(a,"YX") == 0){
    return 39;
  }
  if(strcmp(a,"ZrX") == 0){
    return 40;
  }
  if(strcmp(a,"NbX") == 0){
    return 41;
  }
  if(strcmp(a,"MoX") == 0){
    return 42;
  }
  if(strcmp(a,"TcX") == 0){
    return 43;
  }
  if(strcmp(a,"RuX") == 0){
    return 44;
  }
  if(strcmp(a,"RhX") == 0){
    return 45;
  }
  if(strcmp(a,"PdX") == 0){
    return 46;
  }
  if(strcmp(a,"AgX") == 0){
    return 47;
  }
  if(strcmp(a,"CdX") == 0){
    return 48;
  }
  if(strcmp(a,"InX") == 0){
    return 49;
  }
  if(strcmp(a,"SnX") == 0){
    return 50;
  }
  if(strcmp(a,"SbX") == 0){
    return 51;
  }
  if(strcmp(a,"TeX") == 0){
    return 52;
  }
  if(strcmp(a,"IX") == 0){
    return 53;
  }
  if(strcmp(a,"XeX") == 0){
    return 54;
  }
  if(strcmp(a,"CsX") == 0){
    return 55;
  }
  if(strcmp(a,"BaX") == 0){
    return 56;
  }
  if(strcmp(a,"LaX") == 0){
    return 57;
  }
  if(strcmp(a,"CeX") == 0){
    return 58;
  }
  if(strcmp(a,"PrX") == 0){
    return 59;
  }
  if(strcmp(a,"NdX") == 0){
    return 60;
  }
  if(strcmp(a,"PmX") == 0){
    return 61;
  }
  if(strcmp(a,"SmX") == 0){
    return 62;
  }
  if(strcmp(a,"EuX") == 0){
    return 63;
  }
  if(strcmp(a,"GdX") == 0){
    return 64;
  }
  if(strcmp(a,"TbX") == 0){
    return 65;
  }
  if(strcmp(a,"DyX") == 0){
    return 66;
  }
  if(strcmp(a,"HoX") == 0){
    return 67;
  }
  if(strcmp(a,"ErX") == 0){
    return 68;
  }
  if(strcmp(a,"TmX") == 0){
    return 69;
  }
  if(strcmp(a,"YbX") == 0){
    return 70;
  }
  if(strcmp(a,"LuX") == 0){
    return 71;
  }
  if(strcmp(a,"HfX") == 0){
    return 72;
  }
  if(strcmp(a,"TaX") == 0){
    return 73;
  }
  if(strcmp(a,"WX") == 0){
    return 74;
  }
  if(strcmp(a,"ReX") == 0){
    return 75;
  }
  if(strcmp(a,"OsX") == 0){
    return 76;
  }
  if(strcmp(a,"IrX") == 0){
    return 77;
  }
  if(strcmp(a,"PtX") == 0){
    return 78;
  }
  if(strcmp(a,"AuX") == 0){
    return 79;
  }
  if(strcmp(a,"HgX") == 0){
    return 80;
  }
  if(strcmp(a,"TlX") == 0){
    return 81;
  }
  if(strcmp(a,"PbX") == 0){
    return 82;
  }
  if(strcmp(a,"BiX") == 0){
    return 83;
  }
  if(strcmp(a,"PoX") == 0){
    return 84;
  }
  if(strcmp(a,"AtX") == 0){
    return 85;
  }
  if(strcmp(a,"RnX") == 0){
    return 86;
  }
  if(strcmp(a,"FrX") == 0){
    return 87;
  }
  if(strcmp(a,"RaX") == 0){
    return 88;
  }
  if(strcmp(a,"AcX") == 0){
    return 89;
  }
  if(strcmp(a,"ThX") == 0){
    return 90;
  }
  if(strcmp(a,"PaX") == 0){
    return 91;
  }
  if(strcmp(a,"UX") == 0){
    return 92;
  }
  if(strcmp(a,"NpX") == 0){
    return 93;
  }
  if(strcmp(a,"PuX") == 0){
    return 94;
  }
  if(strcmp(a,"AmX") == 0){
    return 95;
  }
  if(strcmp(a,"CmX") == 0){
    return 96;
  }
  return 1000;
}

////////// << function 7 "Atomicmass" >> ////////// return atomic mass (9410 added)

double Atomicmass(char *a){
  if(strcmp(a,"HX") == 0){
    return 1.008;
  }
  if(strcmp(a,"HeX") == 0){
    return 4.003;
  }
  if(strcmp(a,"LiX") == 0){
    return 6.941;
  }
  if(strcmp(a,"BeX") == 0){
    return 9.012;
  }
  if(strcmp(a,"BX") == 0){
    return 10.81;
  }
  if(strcmp(a,"CX") == 0){
    return 12.01;
  }
  if(strcmp(a,"NX") == 0){
    return 14.01;
  }
  if(strcmp(a,"OX") == 0){
    return 16.00;
  }
  if(strcmp(a,"FX") == 0){
    return 19.00;
  }
  if(strcmp(a,"NeX") == 0){
    return 20.18;
  }
  if(strcmp(a,"NaX") == 0){
    return 22.99;
  }
  if(strcmp(a,"MgX") == 0){
    return 24.31;
  }
  if(strcmp(a,"AlX") == 0){
    return 26.98;
  }
  if(strcmp(a,"SiX") == 0){
    return 28.09;
  }
  if(strcmp(a,"PX") == 0){
    return 30.97;
  }
  if(strcmp(a,"SX") == 0){
    return 32.07;
  }
  if(strcmp(a,"ClX") == 0){
    return 35.45;
  }
  if(strcmp(a,"ArX") == 0){
    return 39.95;
  }
  if(strcmp(a,"KX") == 0){
    return 39.10;
  }
  if(strcmp(a,"CaX") == 0){
    return 40.08;
  }
  if(strcmp(a,"ScX") == 0){
    return 44.96;
  }
  if(strcmp(a,"TiX") == 0){
    return 47.88;
  }
  if(strcmp(a,"VX") == 0){
    return 50.94;
  }
  if(strcmp(a,"CrX") == 0){
    return 52.00;
  }
  if(strcmp(a,"MnX") == 0){
    return 54.94;
  }
  if(strcmp(a,"FeX") == 0){
    return 55.85;
  }
  if(strcmp(a,"CoX") == 0){
    return 58.93;
  }
  if(strcmp(a,"NiX") == 0){
    return 58.69;
  }
  if(strcmp(a,"CuX") == 0){
    return 63.55;
  }
  if(strcmp(a,"ZnX") == 0){
    return 65.39;
  }
  if(strcmp(a,"GaX") == 0){
    return 69.72;
  }
  if(strcmp(a,"GeX") == 0){
    return 72.61;
  }
  if(strcmp(a,"AsX") == 0){
    return 74.92;
  }
  if(strcmp(a,"SeX") == 0){
    return 78.96;
  }
  if(strcmp(a,"BrX") == 0){
    return 79.90;
  }
  if(strcmp(a,"KrX") == 0){
    return 83.80;
  }
  if(strcmp(a,"RbX") == 0){
    return 85.47;
  }
  if(strcmp(a,"SrX") == 0){
    return 87.62;
  }
  if(strcmp(a,"YX") == 0){
    return 88.91;
  }
  if(strcmp(a,"ZrX") == 0){
    return 91.22;
  }
  if(strcmp(a,"NbX") == 0){
    return 92.91;
  }
  if(strcmp(a,"MoX") == 0){
    return 95.94;
  }
  if(strcmp(a,"TcX") == 0){
    return 99.00;
  }
  if(strcmp(a,"RuX") == 0){
    return 101.10;
  }
  if(strcmp(a,"RhX") == 0){
    return 102.90;
  }
  if(strcmp(a,"PdX") == 0){
    return 106.40;
  }
  if(strcmp(a,"AgX") == 0){
    return 107.90;
  }
  if(strcmp(a,"CdX") == 0){
    return 112.40;
  }
  if(strcmp(a,"InX") == 0){
    return 114.80;
  }
  if(strcmp(a,"SnX") == 0){
    return 118.70;
  }
  if(strcmp(a,"SbX") == 0){
    return 121.80;
  }
  if(strcmp(a,"TeX") == 0){
    return 127.60;
  }
  if(strcmp(a,"IX") == 0){
    return 126.90;
  }
  if(strcmp(a,"XeX") == 0){
    return 131.30;
  }
  if(strcmp(a,"CsX") == 0){
    return 132.90;
  }
  if(strcmp(a,"BaX") == 0){
    return 137.30;
  }
  if(strcmp(a,"LaX") == 0){
    return 138.90;
  }
  if(strcmp(a,"CeX") == 0){
    return 140.10;
  }
  if(strcmp(a,"PrX") == 0){
    return 140.90;
  }
  if(strcmp(a,"NdX") == 0){
    return 144.20;
  }
  if(strcmp(a,"PmX") == 0){
    return 145.00;
  }
  if(strcmp(a,"SmX") == 0){
    return 150.40;
  }
  if(strcmp(a,"EuX") == 0){
    return 152.00;
  }
  if(strcmp(a,"GdX") == 0){
    return 157.30;
  }
  if(strcmp(a,"TbX") == 0){
    return 158.90;
  }
  if(strcmp(a,"DyX") == 0){
    return 162.50;
  }
  if(strcmp(a,"HoX") == 0){
    return 164.90;
  }
  if(strcmp(a,"ErX") == 0){
    return 167.30;
  }
  if(strcmp(a,"TmX") == 0){
    return 168.90;
  }
  if(strcmp(a,"YbX") == 0){
    return 173.00;
  }
  if(strcmp(a,"LuX") == 0){
    return 175.00;
  }
  if(strcmp(a,"HfX") == 0){
    return 178.50;
  }
  if(strcmp(a,"TaX") == 0){
    return 180.90;
  }
  if(strcmp(a,"WX") == 0){
    return 183.80;
  }
  if(strcmp(a,"ReX") == 0){
    return 186.20;
  }
  if(strcmp(a,"OsX") == 0){
    return 190.20;
  }
  if(strcmp(a,"IrX") == 0){
    return 192.20;
  }
  if(strcmp(a,"PtX") == 0){
    return 195.10;
  }
  if(strcmp(a,"AuX") == 0){
    return 197.00;
  }
  if(strcmp(a,"HgX") == 0){
    return 200.60;
  }
  if(strcmp(a,"TlX") == 0){
    return 204.40;
  }
  if(strcmp(a,"PbX") == 0){
    return 207.20;
  }
  if(strcmp(a,"BiX") == 0){
    return 209.00;
  }
  if(strcmp(a,"PoX") == 0){
    return 210.00;
  }
  if(strcmp(a,"AtX") == 0){
    return 210.00;
  }
  if(strcmp(a,"RnX") == 0){
    return 222.00;
  }
  if(strcmp(a,"FrX") == 0){
    return 223.00;
  }
  if(strcmp(a,"RaX") == 0){
    return 226.00;
  }
  if(strcmp(a,"AcX") == 0){
    return 227.00;
  }
  if(strcmp(a,"ThX") == 0){
    return 232.00;
  }
  if(strcmp(a,"PaX") == 0){
    return 231.00;
  }
  if(strcmp(a,"UX") == 0){
    return 238.00;
  }
  if(strcmp(a,"NpX") == 0){
    return 237.00;
  }
  if(strcmp(a,"PuX") == 0){
    return 239.00;
  }
  if(strcmp(a,"AmX") == 0){
    return 243.00;
  }
  if(strcmp(a,"CmX") == 0){
    return 247.00;
  }
  if(strcmp(a,"BkX") == 0){
    return 247.00;
  }
  if(strcmp(a,"CfX") == 0){
    return 252.00;
  }
  if(strcmp(a,"EsX") == 0){
    return 252.00;
  }
  if(strcmp(a,"FmX") == 0){
    return 257.00;
  }
  if(strcmp(a,"MdX") == 0){
    return 256.00;
  }
  if(strcmp(a,"NoX") == 0){
    return 259.00;
  }
  if(strcmp(a,"LrX") == 0){
    return 260.00;
  }
  if(strcmp(a,"TVX") == 0 || strcmp(a,"TvX") == 0){
    return 0.000;
  }
  return 0.000;
}
