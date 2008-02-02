function  [n_del,delta,bin_size,n_bin,n_call,bins,pdf]=read1pdf(fid)

% number of deltas computed
n_del=fread(fid,1,'float64');
delta=fread(fid,n_del,'float64');


bin_size=fread(fid,1,'float64');
if (bin_size == 0)
  % this means that we have a bin size for every delta
  bin_size=fread(fid,n_del,'float64');
else
  bin_size=bin_size*ones(n_del);
end 

n_bin=fread(fid,1,'float64');
n_call=fread(fid,1,'float64');

bins=zeros([2*n_bin+1,n_del]);
for i=1:n_del
    bins(:,i)=((-n_bin:n_bin)*bin_size(i))';
end

pdf=fread(fid,[(2*n_bin+1),n_del],'float64');

i=2; % always leave the first one 
while (i<=n_del & n_call>0)
  x=sum(pdf(:,i));
  if abs(x)<1e-9
    % no pdf data for this value of delta, ignore it:
    % disp(sprintf('removing delta =      %i  val=%i',i,delta(i)))
    inew=[1:(i-1),(i+1):n_del];
    delta=delta(inew);
    pdf=pdf(:,inew);
    n_del=n_del-1;
  elseif abs(x-1)>1e-9
    disp(sprintf('warning: pdf sums to: %f',x));
    i=i+1; 
  else
    i=i+1;
  end
end

% return total number of bins:
n_bin=2*n_bin+1;
