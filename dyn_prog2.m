%
% Adam Goldstein
% 2/5/20
%

% Additive dynamics with repeated gamble choice.
% Bet size = k*D, where k=1...K.
% Probability of win is variable p.
% N = #repeated gambles, w0 = initial wealth
% D = minimum bet size (basic grid size of binomial tree)
% K multiplies D for bet size.
% S = put option strike price.  0 means no put option.

global w0 N p eta K n w0_D p_indif k D S base;
global w_array U_array choice_array choice_mat kernel_array eta_d_vec;


function default_global_vars()
global w0 N p eta K D S
  w0 = 4;
  N = 21;
  p = .5498345477415306;
  eta = 1;
  K = 1;
  D = 1;
  S = 0; % default is no put option
endfunction


function display_global_vars()
global w0 N p eta K D S;
  printf("w0=%d D=%d p=%1.8f N=%d eta=%1.2f K=%d S=%1.6f \n", 
          w0, D, p, N, eta, K, S);
endfunction


% Utility of wealth (CRRA), with put option on wealth at strike price S.
% If S=0 it's as if there's no put option since wealth is 
% always positive in that case (absorbing barrier at w=1).
function retval = U(wealth, eta)
global S;
  if (eta==1)
    retval = log(max(wealth,S));
  else
    retval = (max(wealth,S).^(1-eta) - 1)/(1-eta);
  endif
endfunction


function gen_w_array()
global w0 N w_array K D S;
  if ((S==0) && (w0 <= K*D))
    printf("Warning: kernel K at w0 is truncated. \n");
  endif 
  lo = hi = w0;
  for i=1:N+1
    w_array{i} = lo:D:hi;
    hi = hi+K*D;
    if ((S>0) || (lo>K*D))
      lo = lo-K*D;
    endif
  endfor
endfunction


function gen_kernel_array()
global K kernel_array p;
  for k=1:K
    kernel_array{k} = [p zeros(1,2*k-1) 1-p];
  endfor
endfunction


% Dynamic programming using backward recursion
% Generates w_array, choice_array, and kernel_array as a side-effect.
function gen_U_array()
global p N w_array U_array choice_array eta K kernel_array;
  gen_w_array();
  gen_kernel_array();
  U_array{N+1} = U(w_array{N+1},eta);
  %printf("U(N+1): \n");
  %U(w_array{N+1},eta) % Debug
  for n=N:-1:1
    %printf("n=%d \n",n);
    %printf("Initial U(n): \n");
    % debug_U_array_n = U_array{n}
    ind_offset = length(w_array{n+1})-length(w_array{n})-K;
    %for k=1:K
    for k=(K-1):K
      %printf("k=%d \n",k);
      if (k==0)
        Uconv = U_array{n+1};
      else
        Uconv = conv(U_array{n+1},kernel_array{k});
      endif
      ind_lo = (2*k+1) + max(ind_offset-k,0);
      ind_hi = length(w_array{n+1}) - (K-k);
      pad = max(k-ind_offset,0);
      Uconv_adj = [-realmax*ones(1,pad) Uconv(ind_lo:ind_hi)];
      if (k==(K-1))
        accept = ones(1,length(w_array{n}));
        choice_array{n} = zeros(1,length(w_array{n})); % value doesn't matter
        U_array{n} = zeros(1,length(w_array{n})); % value doesn't matter
      else
        accept = Uconv_adj>U_array{n};
      endif
      choice_array{n} = accept.*k + ~accept.*choice_array{n};
      % debug_choice_array_n = choice_array{n}
      U_array{n} = accept.*Uconv_adj + ~accept.*U_array{n};
      % debug_U_array_n = U_array{n}
      % printf("press any key to continue...\n");
      % kbhit();
    endfor
  endfor  
endfunction


% Helper function for calc_eta_d. Computes difference in expected utility
% of gambles k and (k-1) for the given value of p_indif_d. 
function retval = f_p_indif_d(p_indif_d)
global p K k U_array;
  w0_index = K+1;
  p_save = p;
  p = p_indif_d;
  gen_U_array();
  retval = ... 
    p_indif_d*U_array{2}(w0_index+k)+(1-p_indif_d)*U_array{2}(w0_index-k) - ...
   (p_indif_d*U_array{2}(w0_index+k-1)+(1-p_indif_d)*U_array{2}(w0_index-k+1));
  p = p_save;
  %printf("k=%d K=%d p_indif_d=%f retval=%f \n",k,K,p_indif_d,retval);
endfunction


% Helper function for calc_eta_d.  Solves for eta_d that matches
% expected utility of gambles k and (k-1) using previously computed
% win probability p_indif_d.
function retval = f_eta_d(eta_d)
global p_indif_d w0 D k;
  retval = p_indif_d*U(w0+k*D,eta_d)+(1-p_indif_d)*U(w0-k*D,eta_d) -...
           (p_indif_d*U(w0+(k-1)*D,eta_d)+(1-p_indif_d)*U(w0-(k-1)*D,eta_d));
  %printf("eta_d=%f retval=%f k=%d p_indif_d=%f \n",eta_d,retval,k,p_indif_d);
endfunction


% Calculate "derived eta" (eta_d) using current values for w0,D,N and eta. 
% Currently assumes gamble choice is between K and (K-1)l
% eta_d is calculated by first solving for "derived indifference probability"
% p_indif_d, where p=p_indif_d causes accept/reject indifference in period 1.
% 
function eta_d = calc_eta_d()
global N k K p_indif_d;
  if (N==1)
    printf("Warning: N should be greater than 1 for this function. \n");
  endif 
  % Solve for win probabillity p_indif_d that matches expected 
  % derived utility of gamble K to that of (K-1) in period 1.  
  % Calls gen_U_array() repeatedly until matching p_indif_d is found.
  k = K; 
  [p_indif_d, fval, info, output] = fzero(@f_p_indif_d,[0.5,1]);
  % printf("p_indif_d=%1.6f \n", p_indif_d);
  % Solve for eta_d that matches expected utility of gambles k and (k-1) 
  % using previously computed p_indif_d.
  [eta_d, fval, info, output] = fzero(@f_eta_d,[0,1]);
  printf("N=%d p_indif_d=%1.6f eta_d=%1.3f \n",
          N, p_indif_d, eta_d);
endfunction


function plot_eta_d_vec()
global w0 eta D N eta_d_vec p_indif_d S K;
  if (mod(N,2)==0 || N<3)
    printf("Warning: N must be odd and at least equal to 3 \n");
  endif
  eta_d_vec = ones((N-1)/2+1,1);
  eta_d_vec(1)=eta; % by definition
  N_save = N;
  for N=3:2:N_save  % Seems to be problem with even N when S>0
    eta_d_vec((N-1)/2 + 1) = calc_eta_d();
  endfor
  N = N_save;
  plot(1:2:N,eta_d_vec);
  xlabel("Periods");
  ylabel("\\eta_{d}");
  title(sprintf("w0=%1.3f D=%d K=%d eta=%1.2f S=%1.6f\n",w0,D,K,eta,S));
  %title(sprintf("Put option: W_{0}=%1.1f, D=%d, W_{min}=%1.1f \n",w0,D,S));
endfunction


function gen_choice_mat()
global w0 w_array K D N choice_array choice_mat S base;
  hi = w0 + K*D*N;
  if (S==0)
    if (mod(w0,K*D)==0)
      numKD = w0/(K*D) - 1;
    else
      numKD = floor(w0/(K*D));
    endif
    base = 0;
  else
    lo = w0 - K*D*N;
    if (lo<1)
      base = 1 + abs(lo);
    else
      base = 0; % Put is never exercised
    endif
  endif
  choice_mat = (K+1)*ones(base+hi,N); 
  for n=1:N
    for i=1:length(w_array{n})
      choice_mat(w_array{n}(i)+base,n) = choice_array{n}(i);  
    endfor
  endfor
endfunction


function plot_choice_mat()
global choice_mat p N w0 eta K D S;
  if (K==1)
    cmap = [0 0 0;0 1 0;1 1 1];
  elseif (K==2)
    cmap = [0 0 0;0 1 0;0 0 1;1 1 1];
  else
    cmap = [0 0 0; viridis(K);1 1 1];
  endif
  colormap(cmap);
  imagesc(choice_mat+1);
  ax = gca();
  set(ax,"ydir","normal");
  xlabel('Period');
  ylabel('Wealth');
  title(sprintf("w0=%d D=%d p=%1.8f N=%d eta=%1.2f K=%d S=%1.6f \n",
                w0, D, p, N, eta, K, S));
  colorbar;  
endfunction


% Calculate single-period indifference prob for kernel k compared to k-1
% at period n (out of N total).
%
% This function is not called by any others within this file.
function retval = calc_p_indif(k)
global w0 eta D;
    retval = (U(w0-(k-1)*D,eta) - U(w0-k*D,eta)) / ...
             (U(w0+k*D,eta)-U(w0-k*D,eta)-U(w0+(k-1)*D,eta)+U(w0-(k-1)*D,eta));            
endfunction

