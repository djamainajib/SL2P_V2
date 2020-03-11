function [DATA]=outpout_out_of_range_flag_function_SL2P (DATA,box)

   DATA(find(DATA(:,3)<box.Pmin & DATA(:,3)>=box.Pmin-box.Tolerance),5)=1;
   DATA(find(DATA(:,3)>box.Pmax & DATA(:,3)<=box.Pmax+box.Tolerance),6)=1;
   
   DATA(find(DATA(:,3)<box.Pmin-box.Tolerance),7)=1;
   DATA(find(DATA(:,3)>box.Pmax+box.Tolerance),8)=1;
   
   DATA=DATA(:,end-3:end);
end