function Gamma = getChristoffelSym(gu,diff_1_gl,i,k,l)
    Gamma = 0;    
    for m = 1:4
        Gamma = Gamma + 1/2*gu{i,m}.*(diff_1_gl{m,k,l}+diff_1_gl{m,l,k}-diff_1_gl{k,l,m});
    end
end