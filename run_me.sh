for i in m_Bae_50 m_Bea_50 m_ae_50 m_ea_50; do 
    for suf in my orig; do 
	if [ "$suf" == "orig" ]; then
	    ./rcm_${suf}.out input2/${i}.cm 1 > output2/${i}.fm;
	fi
	./rcm_${suf}.out input2/${i}.cm > output2/${i}.res_${suf}.fm;
    done;
done
