%Read MFR Data

for itag = 1:ntags
    for idhvap = 1:ndh     
        for iaccom = 1:naccom
            for iday = 1:ndays
            
            scen2 = [cell2mat(age),'/', cell2mat(tags(itag)),'/dHvap=' num2str(dhvaps(idhvap))];
            fname = ['../Figs/', scen2, '/', 'td.data.alpha_', num2str(accom(iaccom))];
            [day_tmp,hour_tmp,mfr_tmp] = textread(fname, '%s %f %f','headerlines',1);
            
            if strmatch(days(iday),'avg')
                for ihr = 1:24
                    indx1 = ihr;
                    indx2 = ihr+24;
                    indx3 = ihr+48;
                    mfr(itag,idhvap,iaccom,ihr,iday) = ...
                        nanmean([mfr_tmp(indx1),mfr_tmp(indx2),mfr_tmp(indx3)]);
                end
            else
                for jday = 1:ndays
                    indx = strmatch(cell2mat(days(jday)),day_tmp);
                    mfr(itag,idhvap,iaccom,1:24,iday) = mfr_tmp(indx);
                end
            end
            end
        end
    end
end
