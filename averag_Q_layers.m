fid=fopen('Qaveraged_Sarb.txt','w');
top_Lay=[8807 8571 8107 7590 7544 ,7005 6685 6514 6381 6375 ,5415 4538 4469 2882 2162 ,2133 1285 482];
Lay=fliplr(top_Lay);
names_top={'kharaib' 'shuaib' 'nahr-umr' 'mauddud' 'shilaif' 'tuwayil' 'ruwaydha' 'laffan' 'halul' 'fiqa'  'simsima' 'uer-basal-shale' 'umm-er-radhuma' 'rus' 'dammam-basal-shale' 'dammam' 'meocene'};
names=fliplr(names_top);
for i=1:numel(names)
    n=0;sum1=0;sum2=0;sum3=0;
    for j=1:length(Vsc2)
        if(Vsc2(j,1)>=Lay(i) && Vsc2(j,1)<=Lay(i+1))
            sum_vp=sum1+Qtot(j);
            sum_ss=sum2+Qscat(j);
            sum_rho=sum3+(Qtot(j)-Qscat(j));
            n=n+1;
        end;
    end;
    av1=sum1/n;av2=sum2/n;av3=sum3/n;
    fprintf(fid,'\n%f\t%f\t%f\t%s',av1,av2,av3,names{i});
end;