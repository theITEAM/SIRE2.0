void traceinit()                                                       // Initialises the trace plot
{
  long p, c, cl, f, z;

  cout << "1|";

  if(snpfl == 1){
    cout << "SNP|a_g|SNP effect for susceptibility|";
    cout << "SNP|a_f|SNP effect for infectivity|";
    if(mod == SIR) cout << "SNP|a_r|SNP effect for recoverability|";
    if(domon == 1){
      cout << "SNP|Δ_g|Dominance factor for susceptibility|";
      cout << "SNP|Δ_f|Dominance factor for infectivity|";
      if(mod == SIR) cout << "SNP|Δ_r|Dominance factor for recoverability|";
    }
  }

  for(f = 0; f < nfi; f++){
    cout << "Fix. Eff.|" << finame[0][f] << "|Fixed effect " << finame[0][f] << " for suceptibility|";
    cout << "Fix. Eff.|" << finame[1][f] << "|Fixed effect " << finame[1][f] << " for infectivity|";
    if(mod == SIR) cout << "Fix. Eff.|" << finame[2][f] << "|Fixed effect " << finame[2][f] << " for recoverability|";
  }

  cout << "Epi.|β|Contact rate|";
  if(mod == SIR) cout << "Epi.|γ|Recovery rate|";
  if(mod == SIR) cout << "Epi.|k|Shape parameter|";

  if(envon == 1){
    cout << "Covar.|Ψ_gg|Covariance matrix Ψ_gg for residule|";
    cout << "Covar.|Ψ_gf|Covariance matrix Ψ_gf for residule|";
    if(mod == SIR) cout << "Covar.|Ψ_gr|Covariance matrix Ψ_gr for residule|";
    cout << "Covar.|Ψ_ff|Covariance matrix Ψ_ff for residule|";
    if(mod == SIR) cout << "Covar.|Ψ_fr|Covariance matrix Ψ_fr for residule|";
    if(mod == SIR) cout << "Covar.|Ψ_rr|Covariance matrix Ψ_rr for residule|";
  }

  if(randon == 1){
    cout << "Covar.|Ω_gg|Covariance matrix Ω_gg for random effect|";
    cout << "Covar.|Ω_gf|Covariance matrix Ω_gf for random effect|";
    if(mod == SIR) cout << "Covar.|Ω_gr|Covariance matrix Ω_gr for random effect|";
    cout << "Covar.|Ω_ff|Covariance matrix Ω_ff for random effect|";
    if(mod == SIR) cout << "Covar.|Ω_fr|Covariance matrix Ω_fr for random effect|";
    if(mod == SIR) cout << "Covar.|Ω_rr|Covariance matrix Ω_rr for random effect|";
  }

  if(envon == 1 && randon == 1){
    cout << "Herit.|h2_g|Heritability of susceptibility|";
    cout << "Herit.|h2_f|Heritability of infectivity|";
    if(mod == SIR) cout << "Herit.|h2_r|Heritability of recoverability|";
  }


  if(geffon == 1 && Z > 1){
    cout << "Gr. Eff.|σ_G|Standard deviation in group effect" << "|";
    for(z = 0; z < Z; z++) cout << "Gr. Eff.|G" << z << "|Group effect for trial " << z << "|";
  }


  cout << "Misc.|L_inf|Infection likelihood|";
  if(mod == SIR) cout << "Misc.|L_rec|Recovery likelihood|";
  cout << "Misc.|L_DT|Likelihood of diagnostic tests|";
  cout << "Misc.|L_e|Environmental likelihood|";
  cout << "Misc.|L_G|Group effect likelihood|";
  cout << "Misc.|Pr|Prior|";


	for(int i = 0; i < predacc.size(); i++){
		if(q_g_bv.size() > 0){
			cout << "Pred. Ac.|"+predacc[i].name+": Sus.|"+predacc[i].name+" for sus.|";
		}
		
		if(q_f_bv.size() > 0){
			cout << "Pred. Ac.|"+predacc[i].name+": Inf.|"+predacc[i].name+" for inf.|";
		}
		
		if(q_r_bv.size() > 0){
			cout << "Pred. Ac.|"+predacc[i].name+": Rec.|"+predacc[i].name+" for rec.|";
		}
	}
    
  cout << "\n";
  cout.flush();
}

void traceplot()                                                      // Outputs varaible values
{
  long f, z;

  cout <<"0|";
  if(snpfl == 1){
    cout << a_g << "|" << a_f << "|"; if(mod == SIR) cout << a_r << "|";
    if(domon == 1){
      cout << delta_g << "|" << delta_f << "|"; if(mod == SIR) cout << delta_r << "|";
    }
  }

  for(f = 0; f < nfi; f++){ cout << fi_g[f] << "|" << fi_f[f] << "|"; if(mod == SIR) cout << fi_r[f] << "|";}

  cout << beta << "|"; if(mod == SIR) cout << gama << "|" << kshape << "|";

  if(envon == 1){
    if(mod == SIR) cout << vare_gg << "|" << vare_gf << "|" << vare_gr << "|" << vare_ff << "|" << vare_fr << "|" << vare_rr << "|";
    else cout << vare_gg << "|" << vare_gf << "|" << vare_ff << "|";
  }

  if(randon == 1){
    if(mod == SIR) cout << vara_gg << "|" << vara_gf << "|" << vara_gr << "|" << vara_ff << "|" << vara_fr << "|" << vara_rr << "|";
    else cout << vara_gg << "|" << vara_gf << "|" << vara_ff << "|";
  }

  if(envon == 1 && randon == 1){
    cout << vara_gg/(vara_gg+vare_gg) << "|";
    cout << vara_ff/(vara_ff+vare_ff) << "|";
    if(mod == SIR) cout << vara_rr/(vara_rr+vare_rr) << "|";
  }


  if(geffon == 1 && Z > 1){
    cout << siggeff_g << "|";
    for(z = 0; z < Z; z++) cout << geff_g[z] << "|";
  }

  cout << Li << "|";
  if(mod == SIR) cout << Lrec << "|";

  cout << Ldtest << "|" << Lie << "|" << Ligeff_g << "|" << Pri << "|";
  
	if(nqsum < 2){
		for(int pa = 0; pa < predacc.size(); pa++){
			if(q_g_bv.size() > 0) cout << "0|";
			if(q_f_bv.size() > 0) cout << "0|";
			if(mod == SIR && q_r_bv.size() > 0) cout << "0|";
		}
	}
	else{
		vector <double> qgav = q_g_sum, qfav = q_f_sum, qrav = q_r_sum;
	
		for(int i = 0; i < N; i++){
			qgav[i] /= nqsum;
			qfav[i] /= nqsum;
			if(mod == SIR) qrav[i] /= nqsum;
		}

		for(int pa = 0; pa < predacc.size(); pa++){
			if(q_g_bv.size() > 0) cout << correlation(q_g_bv,qgav,predacc[pa].ind) << "|";
			if(q_f_bv.size() > 0) cout << correlation(q_f_bv,qfav,predacc[pa].ind) << "|";
			if(mod == SIR && q_r_bv.size() > 0) cout << correlation(q_r_bv,qrav,predacc[pa].ind) << "|";
		}
	}
	
  cout << "\n";
  cout.flush();
}

double correlation(const vector <double> &val1, const vector <double> &val2, const vector <long> &ind)
{
	double av1 = 0, av2 = 0, av11 = 0, av22 = 0, av12 = 0;
  
	int jmax = ind.size();
	for(int j = 0; j < jmax; j++){
		int i = ind[j];
	
		av1 += val1[i]; av2 += val2[i]; 
		av11 += val1[i]*val1[i]; av22 += val2[i]*val2[i];
		av12 += val1[i]*val2[i];  
  }
	av1 /= jmax; av2 /= jmax; av11 /= jmax; av22 /= jmax; av12 /= jmax;

  return (av12 - av1*av2)/(sqrt(av11 - av1*av1)*sqrt(av22 - av2*av2));
}

void eventplot()                                            // Outputs infection and recovery events
{
  long i, z, k;
  vector <long> list;
  vector <double> listt;

  cout << "5|";

  cout << N << "|";
  for(i = 0; i < N; i++){
    z = indtrial[i];
    
    list.clear(); listt.clear();
    if(z >= 0){
      list.push_back(0); listt.push_back(trialtmin[z]);

      if(It[i] != -big){
	list.push_back(1); listt.push_back(It[i]);
	if(Rt[i] < trialtmax[z]){ list.push_back(2); listt.push_back(Rt[i]);}
      }
      list.push_back(-1); listt.push_back(trialtmax[z]);

      for(k = 0; k < list.size()-1; k++) if(listt[k] > listt[k+1]) emsg("Order of events not right");
    }
    
    cout << list.size() << "|"; for(k = 0; k < list.size(); k++) cout << list[k] << "|" << listt[k] << "|";
  }
  cout << "\n";
}

void diagnostic()                                                          // Outputs MCMC diagnostics
{
  long i, fi, ty;
  double timeto = timetot + clock(), fav, fmin, fmax, f, nfav;

  cout << "6|";
  cout <<  "Acceptance Ratios|";
  for(i = 0; i < 5; i++) cout << nac[i]/ntr[i] << " " ; cout << " Parameters|";
  if(mod == SIR){ for(i = 0; i < 4; i++) cout << nac_rec[i]/ntr_rec[i] << " " ; cout << " Rec Param|";}

  if(envon == 1){
    cout << nac_e_g/ntr_e_g << " e_g   ";
    cout << nac_e_f/ntr_e_f << " e_f    ";
    if(mod == SIR) cout << nac_e_r/ntr_e_r << " e_r";
    cout << "|";

    for(ty = 0; ty < 9; ty++)  cout << nac_PBP[ty]/ntr_PBP[ty] << " fa:" <<  nfa_PBP[ty]/ntr_PBP[ty] << " " <<jumpPBP[ty] << " PBP   "<< ty << "|";

    cout << nac_vare_gg/ntr_vare_gg << " vare_gg   ";
    cout << nac_vare_gf/ntr_vare_gf << " vare_gf   ";
    if(mod == SIR) cout << nac_vare_gr/ntr_vare_gr << " vare_gr   ";
    cout << nac_vare_ff/ntr_vare_ff << " vare_ff   ";
    if(mod == SIR) cout << nac_vare_fr/ntr_vare_fr << " vare_fr   ";
    if(mod == SIR) cout << nac_vare_rr/ntr_vare_rr << " vare_rr   ";
    cout << "|";
  }

  if(randon == 1){
    cout << nac_q_g/ntr_q_g << " q_g   ";
    cout << nac_q_f/ntr_q_f << " q_f    ";
    if(mod == SIR) cout << nac_q_r/ntr_q_r << " q_r    ";
    cout << "|";
    cout << nac_vara_gg/ntr_vara_gg << " vara_gg   ";
    cout << nac_vara_gf/ntr_vara_gf << " vara_gf   ";
    if(mod == SIR)cout << nac_vara_gr/ntr_vara_gr << " vara_gr   ";
    cout << nac_vara_ff/ntr_vara_ff << " vara_ff   ";
    if(mod == SIR) cout << nac_vara_fr/ntr_vara_fr << " vara_fr   ";
    if(mod == SIR) cout << nac_vara_rr/ntr_vara_rr << " vara_rr";
    cout << "|";
  }

  if(geffon == 1){
    cout << nac_geff_g/ntr_geff_g << " geff_g   ";
    cout << nac_siggeff_g/ntr_siggeff_g << " siggeff_g  |";
  }

  for(fi = 0; fi < nfi; fi++){
    cout << fi << " " << nac_fi_g[fi]/ntr_fi_g[fi] << " fi_g   ";
    cout << fi << " " << nac_fi_f[fi]/ntr_fi_f[fi] << " fi_f   ";
    if(mod == SIR) cout << fi << " " << nac_fi_r[fi]/ntr_fi_r[fi] << " fi_r";
    cout << "|";
  }

  long j;

  fmin = 1; fmax = 0; fav = 0; nfav = 0;
  for(j = 0; j < N; j++){ if(ntr_I[j] > 0){ 
		f = nac_I[j]/ntr_I[j]; if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f; nfav++;}
		
		if(f == 0) cout << "zero:" << indid[j]; 
	}
  cout << "shiftI ac:" << fav/nfav << ", " << fmin << "-" << fmax << "|";

  if(mod == SIR){
    fmin = 1; fmax = 0; fav = 0; nfav = 0;
    for(j = 0; j < N; j++){ if(ntr_R[j] > 0){ f = nac_R[j]/ntr_R[j]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f; nfav++;}}
    cout << "shiftR ac:" << fav/nfav << ", " << fmin << "-" << fmax << "|";
  }

  fmin = 1; fmax = 0; fav = 0; nfav = 0;
  for(j = 0; j < N; j++){ if(Iposst[j] == unknownI && ntr_addinf[j] > 0){ f = nac_addinf[j]/ntr_addinf[j]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f; nfav++;}}
  cout << "addinf ac:" << fav/nfav << ", " << fmin << "-" << fmax << "|";

  fmin = 1; fmax = 0; fav = 0; nfav = 0;
  for(j = 0; j < N; j++){ if(Iposst[j] == unknownI && ntr_reminf[j] > 0){ f = nac_reminf[j]/ntr_reminf[j]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f; nfav++;}}
  cout << "reminf ac:" << fav/nfav << ", " << fmin << "-" << fmax << "|";

  cout << "Times:|";
  cout << long(100*time_param/timeto) << " "  << long(100*time_paraminit/timeto) << " time param|";
  cout << long(100*time_rec/timeto) << " "  << long(100*time_recinit/timeto) << " time rec|";

  if(envon == 1){
    cout << long(100*time_e_g/timeto) << " " << long(100*time_e_ginit/timeto) << " time e_g|";
    cout << long(100*time_e_f/timeto) << " " << long(100*time_e_finit/timeto) << " time e_f|";
    if(mod == SIR) cout << long(100*time_e_r/timeto) << " " << long(100*time_e_rinit/timeto) << " time e_r|";

    cout << long(100*time_eq/timeto) << " time eq|";
    cout << long(100*time_PBPe_g/timeto) << " time PBPe_g  ";
    cout << long(100*time_PBPe_f/timeto) << " time PBPe_f   ";
    if(mod == SIR) cout << long(100*time_PBPe_r/timeto) << " time PBPe_r";
    cout << "|";
  }

  if(randon == 1){
    cout << long(100*time_q_g/timeto) << " " << long(100*time_q_ginit/timeto) << " time q_g|";
    cout << long(100*time_q_f/timeto) << " " << long(100*time_q_finit/timeto) << " time q_f|";
    if(mod == SIR) cout << long(100*time_q_r/timeto) << " " << long(100*time_q_rinit/timeto) << " time q_r|";
    cout <<  long(100*time_vara/timeto) << " time vara|";
  }

  if(geffon == 1){
    cout << long(100* time_geff_g /timeto) << " time greff_g|";
  }

  cout << long(100*time_fi/timeto) << " time fi|";
  cout << long(100*time_events/timeto) << " time events|";
  cout << long(100*time_addreminf/timeto) << " time addreminf|";
  cout << long(100*time_init2/timeto) << " time init2|";

  cout << long(100*time_check/timeto) << " time check|";
  cout << timeto/(60.0*CLOCKS_PER_SEC) << "Total time|";

  cout << "\n";
  cout.flush();
}

struct ParamterStore{    // Stores parameter samples
	string name;
	vector <double> sample;
};

long vout;
vector <ParamterStore> pstore;

void store_sample()                                                    // Stores a sample for the parameter values
{
	long f, z, i;
	
	vout = 0;
	
  if(snpfl == 1){
		addsamp(a_g,"a_g"); addsamp(a_f,"a_f"); if(mod == SIR) addsamp(a_r,"a_r");
    if(domon == 1){
			addsamp(delta_g,"Δ_g"); addsamp(delta_f,"Δ_f"); if(mod == SIR) addsamp(delta_r,"Δ_r"); 
    }
  }

  for(f = 0; f < nfi; f++){
		stringstream ss; ss << f; string st = ss.str();
		addsamp(fi_g[f],"b_g_"+st); addsamp(fi_f[f],"b_f_"+st); if(mod == SIR) addsamp(fi_r[f],"b_r_"+st); 
	}

	addsamp(beta,"β"); if(mod == SIR){ addsamp(gama,"γ"); addsamp(kshape,"k");}
  
  if(envon == 1){
    if(mod == SIR){
			addsamp(vare_gg,"Ψ_gg"); addsamp(vare_gf,"Ψ_gf"); addsamp(vare_gr,"Ψ_gr"); 
			addsamp(vare_ff,"Ψ_ff"); addsamp(vare_fr,"Ψ_fr"); addsamp(vare_ff,"Ψ_rr");
		}
    else{
			addsamp(vare_gg,"Ψ_gg"); addsamp(vare_gf,"Ψ_gf"); addsamp(vare_gr,"Ψ_ff"); 
		}
  }
	
	if(randon == 1){
    if(mod == SIR){
			addsamp(vara_gg,"Ω_gg"); addsamp(vara_gf,"Ω_gf"); addsamp(vara_gr,"Ω_gr"); 
			addsamp(vara_ff,"Ω_ff"); addsamp(vara_fr,"Ω_fr"); addsamp(vara_rr,"Ω_rr");
		}
    else{
			addsamp(vara_gg,"Ω_gg"); addsamp(vara_gf,"Ω_gf"); addsamp(vara_gr,"Ω_ff"); 
		}
  }

  if(envon == 1 && randon == 1){
		addsamp(vara_gg/(vara_gg+vare_gg),"h2_gg"); addsamp(vara_ff/(vara_ff+vare_ff),"h2_ff"); 
		if(mod == SIR) addsamp(vara_rr/(vara_rr+vare_rr),"h2_rr"); 
  }

  if(geffon == 1 && Z > 1){
		addsamp(siggeff_g,"σ_G");
		/*
		for(z = 0; z < Z; z++){
			stringstream ss; ss << "G_" << z;
			addsamp(geff_g[z],ss.str());
		}			 
		*/
  }
}

void addsamp(double value, string name)
{
	if(pstore.size() == vout){
		pstore.push_back(ParamterStore ());
		pstore[vout].name = name;
	}
	pstore[vout].sample.push_back(value);
	vout++;
}

struct Statistics{                                          // Stores statistical information
	string mean;                                        // The mean
	string CImin, CImax;                                // The minimum and maximum of the 90% credible interval
	string ESS;                                         // The estimated sample size
};

string to_string(double st)
{
	stringstream ss; ss << st;
	return ss.str();
}

Statistics get_statastic(const vector <double> &vec)                      
{
	double sum, sum2, var, f, sd, a, cor;
	long i, n, d;
	Statistics stat;
	
	n = vec.size();
	if(n == 0){
		stat.mean = "---"; stat.CImin = "---"; stat.CImax = "---"; stat.ESS = "---"; 
	}
	else{
		sum = 0.0, sum2 = 0.0; for(i = 0; i < vec.size(); i++){ sum += vec[i]; sum2 += vec[i]*vec[i];}
		sum /= n; sum2 /= n;
		stat.mean = to_string(sum); 
		
		vector <double> vec2 = vec;
		sort(vec2.begin(),vec2.end());
	
		if(n >= 2){
			i = (unsigned int)((n-1)*0.025); f = (n-1)*0.025 - i;
			stat.CImin = to_string(vec2[i]*(1-f) + vec2[i+1]*f);
				
			i = (unsigned int)((n-1)*0.975); f = (n-1)*0.975 - i;
			stat.CImax = to_string(vec2[i]*(1-f) + vec2[i+1]*f);
		}
		else{
			stat.CImin = to_string(vec2[0]);
			stat.CImax = to_string(vec2[0]);
		}

		vec2 = vec;
		var = sum2 - sum*sum;
		if(var <= 0.0000000001 || n <= 2)  stat.ESS = "---";
		else{	
			sd = sqrt(var);
			for(i = 0; i < vec2.size(); i++) vec2[i] = (vec2[i]-sum)/sd;
				
			sum = 1.0;
			for(d = 1; d < n/2; d++){             // calculates the effective sample size
				a = 0.0; for(i = 0; i < n-d; i++) a += vec2[i]*vec2[i+d]; 
				cor = a/(n-d);
				if(cor < 0) break;
				sum += 2*cor;			
			}
			stat.ESS = to_string(int(n/sum));
		}
	}
	
	return stat;
}

void output_statistics()                                               // Outputs statistics
{
	long v, i, pa;
	
	cout << "Name\tMean\t95% CI (min - max)\tESS\n";
	for(v = 0; v < pstore.size(); v++){
		Statistics stat = get_statastic(pstore[v].sample);
		cout << pstore[v].name << "\t" << stat.mean << "\t(" << stat.CImin << " - " << stat.CImax << ")\t" << stat.ESS << "\n";
	}
	
	for(i = 0; i < N; i++){ q_g_sum[i] /= nqsum; q_f_sum[i] /= nqsum; q_r_sum[i] /= nqsum;}
	
	cout << "\n";
	for(pa = 0; pa < predacc.size(); pa++){
		if(q_g_bv.size() != 0){
			cout <<  predacc[pa].name << " prediction accuracy for susceptibility = " << correlation(q_g_sum,q_g_bv,predacc[pa].ind) << "\n";
		}
		
		if(q_f_bv.size() != 0){
			cout <<  predacc[pa].name << " prediction accuracy for infectivity = " << correlation(q_f_sum,q_f_bv,predacc[pa].ind) << "\n";
		}
		
		if(q_r_bv.size() != 0){
			cout <<  predacc[pa].name << " prediction accuracy for recoverability = " << correlation(q_r_sum,q_r_bv,predacc[pa].ind) << "\n";
		}
	}
}
