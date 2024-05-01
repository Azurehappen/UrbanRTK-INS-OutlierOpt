function SSR = parser_CNES(CNES_path)


SSRfile = fopen(CNES_path);
obt_i = 0;
clk_i = 0;
cdb_i = 0;
phb_i = 0;
type = 0;
SSR.GPS=[];
SSR.GLO=[];
SSR.GAL=[];
SSR.BDS=[];
while ~feof(SSRfile)
   line = fgetl(SSRfile);
   if strcmp(line(1),'>') % new period
       linesplit = strsplit(line);
       switch linesplit{2}
           case 'ORBIT'
               obt_i = obt_i+1;
               [~,~,SSR.orbit_iTOW(:,obt_i)] = date2gnsst(str2double(linesplit(3:8)));%  GPS seconds
               
               type = 1;
           case 'CLOCK'
               clk_i = clk_i+1;
%                IGS.clk_tr(:,clk_i) = str2double(linesplit(3:8))'; % [year;month;date;hour;minute;second]
               [~,~,SSR.clk_iTOW(:,clk_i)] = date2gnsst(str2double(linesplit(3:8)));%  GPS seconds
               
               type = 2;
           case 'CODE_BIAS'
               cdb_i = cdb_i+1;
               [~,~,SSR.cdb_iTOW(:,cdb_i)] = date2gnsst(str2double(linesplit(3:8)));%  GPS seconds
               type = 3;
           case 'PHASE_BIAS'
               phb_i = phb_i+1;
               [~,~,SSR.phb_iTOW(:,cdb_i)] = date2gnsst(str2double(linesplit(3:8)));%  GPS seconds
               type = 4;
       end
   else
       linesp = [line(1),' ',line(2:end)];
       linesplit = strsplit(linesp);
       i = str2double(linesplit(2));
       switch type
           case 1
               switch linesplit{1}
                   case 'G'
                       SSR.GPS.orbit_IDOE(i,obt_i) = str2double(linesplit(3));
                       SSR.GPS.orbit_x(i,obt_i) = str2double(linesplit(4));
                       SSR.GPS.orbit_y(i,obt_i) = str2double(linesplit(5));
                       SSR.GPS.orbit_z(i,obt_i) = str2double(linesplit(6));
                       SSR.GPS.orbit_xv(i,obt_i) = str2double(linesplit(7));
                       SSR.GPS.orbit_yv(i,obt_i) = str2double(linesplit(8));
                       SSR.GPS.orbit_zv(i,obt_i) = str2double(linesplit(9));
                   case 'R'
                       SSR.GLO.orbit_IDOE(i,obt_i) = str2double(linesplit(3));
                       SSR.GLO.orbit_x(i,obt_i) = str2double(linesplit(4));
                       SSR.GLO.orbit_y(i,obt_i) = str2double(linesplit(5));
                       SSR.GLO.orbit_z(i,obt_i) = str2double(linesplit(6));
                       SSR.GLO.orbit_xv(i,obt_i) = str2double(linesplit(7));
                       SSR.GLO.orbit_yv(i,obt_i) = str2double(linesplit(8));
                       SSR.GLO.orbit_zv(i,obt_i) = str2double(linesplit(9));
                   case 'E'
                       SSR.GAL.orbit_IDOE(i,obt_i) = str2double(linesplit(3));
                       SSR.GAL.orbit_x(i,obt_i) = str2double(linesplit(4));
                       SSR.GAL.orbit_y(i,obt_i) = str2double(linesplit(5));
                       SSR.GAL.orbit_z(i,obt_i) = str2double(linesplit(6));
                       SSR.GAL.orbit_xv(i,obt_i) = str2double(linesplit(7));
                       SSR.GAL.orbit_yv(i,obt_i) = str2double(linesplit(8));
                       SSR.GAL.orbit_zv(i,obt_i) = str2double(linesplit(9));
                   case 'C'
                       SSR.BDS.orbit_IDOE(i,obt_i) = str2double(linesplit(3));
                       SSR.BDS.orbit_x(i,obt_i) = str2double(linesplit(4));
                       SSR.BDS.orbit_y(i,obt_i) = str2double(linesplit(5));
                       SSR.BDS.orbit_z(i,obt_i) = str2double(linesplit(6));
                       SSR.BDS.orbit_xv(i,obt_i) = str2double(linesplit(7));
                       SSR.BDS.orbit_yv(i,obt_i) = str2double(linesplit(8));
                       SSR.BDS.orbit_zv(i,obt_i) = str2double(linesplit(9));
               end
           case 2 
               switch linesplit{1}
                   case 'G'
                       SSR.GPS.clk_IDOC(i,clk_i) = str2double(linesplit(3));
                       SSR.GPS.clk_corr(i,clk_i) = str2double(linesplit(4));
                       SSR.GPS.clk_vel(i,clk_i) = str2double(linesplit(5));
                       SSR.GPS.clk_acc(i,clk_i) = str2double(linesplit(6));
                   case 'R'
                       SSR.GLO.clk_IDOC(i,clk_i) = str2double(linesplit(3));
                       SSR.GLO.clk_corr(i,clk_i) = str2double(linesplit(4));
                       SSR.GLO.clk_vel(i,clk_i) = str2double(linesplit(5));
                       SSR.GLO.clk_acc(i,clk_i) = str2double(linesplit(6));
                   case 'E'
                       SSR.GAL.clk_IDOC(i,clk_i) = str2double(linesplit(3));
                       SSR.GAL.clk_corr(i,clk_i) = str2double(linesplit(4));
                       SSR.GAL.clk_vel(i,clk_i) = str2double(linesplit(5));
                       SSR.GAL.clk_acc(i,clk_i) = str2double(linesplit(6));
                   case 'C'
                       SSR.BDS.clk_IDOC(i,clk_i) = str2double(linesplit(3));
                       SSR.BDS.clk_corr(i,clk_i) = str2double(linesplit(4));
                       SSR.BDS.clk_vel(i,clk_i) = str2double(linesplit(5));
                       SSR.BDS.clk_acc(i,clk_i) = str2double(linesplit(6));
               end
           case 3
               switch linesplit{1}
                   case 'G'
                       for k = 3:length(linesplit)
                           if strcmp(linesplit{k},'1C')
                               SSR.GPS.code_bias_C1C(i,cdb_i) = str2double(linesplit(k+1));
                           end
                           if strcmp(linesplit{k},'2L')
                               SSR.GPS.code_bias_C2L(i,cdb_i) = str2double(linesplit(k+1));
                           end
                       end
                   case 'R'
                       for k = 3:length(linesplit)
                           if strcmp(linesplit{k},'1C')
                               SSR.GLO.code_bias_C1C(i,cdb_i) = str2double(linesplit(k+1));
                           end
                           if strcmp(linesplit{k},'2C')
                               SSR.GLO.code_bias_C2C(i,cdb_i) = str2double(linesplit(k+1));
                           end
                       end
                   case 'E'
                       for k = 3:length(linesplit)
                           if strcmp(linesplit{k},'1X')
                               SSR.GAL.code_bias_C1C(i,cdb_i) = str2double(linesplit(k+1));
                           end
                           if strcmp(linesplit{k},'7X')
                               SSR.GAL.code_bias_C2C(i,cdb_i) = str2double(linesplit(k+1));
                           end
                       end
                   case 'C'
                       for k = 3:length(linesplit)
                           if strcmp(linesplit{k},'2I')
                               SSR.BDS.code_bias_C2I(i,cdb_i) = str2double(linesplit(k+1));
                           end
                           if strcmp(linesplit{k},'7I')
                               SSR.BDS.code_bias_C7I(i,cdb_i) = str2double(linesplit(k+1));
                           end
                       end   
               end
           case 4
               
       end
              
   end
end    
fclose(SSRfile);
end