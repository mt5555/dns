function fid=endianopen(fname,access)
%
%  open a file, read first 8 bytes as a float, check
%  for an integer between -32000 and 32000.
%  if not, try again with endianness swapped
%
fid=fopen(fname,access,'l');
if (fid<0) return; end;

test=fread(fid,1,'float64');
if (floor(test)==test & test>-32000 & test < 32000) 
   fseek(fid,0,-1);
   return; 
end

% try big endian
fclose(fid);
fid=fopen(fname,access,'b');
test=fread(fid,1,'float64');
if (floor(test)==test & test>-32000 & test < 32000) 
   fseek(fid,0,-1);
   return; 
end
return





