%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Title                : PICASSO allows ultra-multiplexed fluorescence imaging of spatially overlapping proteins without reference spectra measurements 
%% Authors              : Junyoung Seo, Yeonbo Sim, Jeewon Kim, Hyunwoo Kim, In Cho, Hoyeon Nam, Young-Gyu Yoon and Jae-Byum Chang
%% Authors' Affiliation : Korea Advanced Institute of Science and Technology
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [imgDemixed, unmixing_log] = PICASSO(input_image, qN, maxIter, step_size, monitor)


imgChannel1 = input_image(:,:,1);
imgChannel2 = input_image(:,:,2);
imgChannel3 = input_image(:,:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binF = 16;
qN = 100;
imgResolution = size(imgChannel1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization
u1_img = double(imgChannel1);
u2_img = double(imgChannel2);
u3_img = double(imgChannel3);
v1_img = u1_img;
v2_img = u2_img;
v3_img = u3_img;

%% initial unmixing matrix
UMM = eye(3);

%% initialize log variables
unmixing_log.mi12 = zeros(1, maxIter);
unmixing_log.mi13 = zeros(1, maxIter);
unmixing_log.mi22 = zeros(1, maxIter);
unmixing_log.nr1 = zeros(1, maxIter);
unmixing_log.nr2 = zeros(1, maxIter);
unmixing_log.nr3 = zeros(1, maxIter);
unmixing_log.max_alpha = zeros(1, maxIter);
unmixing_log.alpha_matrix = zeros(3, 3, maxIter);

if monitor
   h3 = figure(3); set(h3, 'Name', 'monitoring iterations...'); 
end

for it=1:maxIter,
   
    disp([num2str(it) ' | ' num2str(maxIter) 'th iteration is running...']);
        
    alpha_matrix = zeros(3,3); 

    %% part 1    
    alpha1 = estimate_alpha(v1_img, v2_img, binF, qN);
    alpha2 = estimate_alpha(v1_img, v3_img, binF, qN);
    alpha_matrix(1,2) =  - alpha1;
    alpha_matrix(1,3) =  - alpha2;
    
    %% part 2
    alpha1 = estimate_alpha(v2_img, v1_img, binF, qN);
    alpha2 = estimate_alpha(v2_img, v3_img, binF, qN);
    alpha_matrix(2,1) =  - alpha1;
    alpha_matrix(2,3) =  - alpha2;
        
    %% part 3  
    alpha1 = estimate_alpha(v3_img, v1_img, binF, qN);
    alpha2 = estimate_alpha(v3_img, v2_img, binF, qN);
    alpha_matrix(3,1) = - alpha1;
    alpha_matrix(3,2) = - alpha2;

    
    %% construct incremental unmixing matrix   
    cUMM = eye(3) + step_size*alpha_matrix;
    UMM = cUMM*UMM;

    
    %% simuletaneous update    
    v1_img = cUMM(1,1)*v1_img + cUMM(1,2)*v2_img + cUMM(1,3)*v3_img;
    v2_img = cUMM(2,1)*v1_img + cUMM(2,2)*v2_img + cUMM(2,3)*v3_img;
    v3_img = cUMM(3,1)*v1_img + cUMM(3,2)*v2_img + cUMM(3,3)*v3_img;   
   
                
    %% construct current image
    monitor_image(:,:,1) = v1_img;
    monitor_image(:,:,2) = v2_img;
    monitor_image(:,:,3) = v3_img;
    monitor_image = perChNorm(monitor_image);


    %% show current image
    if monitor
        if mod(it,5)==1
            figure(3); imshow(monitor_image,'border','tight'); drawnow;
        end                       
    end
    
    if mod(it,50)==1
        v1_img = max(v1_img,0);
        v2_img = max(v2_img,0);
        v3_img = max(v3_img,0); 
    end
    

    unmixing_log.alpha_matrix(:,:,it) = cUMM;     
    unmixing_log.max_alpha(it)  = max(abs(alpha_matrix(:)));
    unmixing_log.nr1(it) = checkNegativity(monitor_image(:,:,1) );
    unmixing_log.nr2(it) = checkNegativity(monitor_image(:,:,2) );
    unmixing_log.nr3(it) = checkNegativity(monitor_image(:,:,3) );
    
        
    
end

unmixing_log.UMM = UMM;
unmixing_log.it = it;

%% enforce positivity to the final result
imgChannel1orth = max(v1_img, 0);
imgChannel2orth = max(v2_img, 0);
imgChannel3orth = max(v3_img, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Show & save result %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgDemixed = zeros( imgResolution(1), imgResolution(2), 3, 'single');
imgDemixed(:,:,1) = imgChannel1orth;
imgDemixed(:,:,2) = imgChannel2orth;
imgDemixed(:,:,3) = imgChannel3orth;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha = estimate_alpha(u1_img, u2_img, binF, qN)

u1_img_bin = pixelBin(u1_img, binF);
u2_img_bin = pixelBin(u2_img, binF);

u1bin = u1_img_bin(:);
u2bin = u2_img_bin(:);

MI_output = @(lambda) calcMI_v3( u1bin - lambda*u2bin, u2bin, qN);

alpha0 = dot(u1bin,u2bin)/sqrt(dot(u1bin,u1bin)*dot(u2bin,u2bin));
alpha = fminsearch(MI_output, 0.1*alpha0);



alpha = min(alpha, 0.5);
alpha = max(alpha, -0.5);



v1bin = u1bin - alpha*u2bin;
negative_ratio = checkNegativity(v1bin);
if negative_ratio>0.9e-3;
   alpha = alpha/10;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function negative_ratio = checkNegativity(img_in)

img_L1 = sum(abs(img_in(:)));
img_negative = (img_in).*(img_in<0);
img_negative_L1 = sum(abs(img_negative(:)));

negative_ratio = img_negative_L1/img_L1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imgOut = pixelBin(imgIn, binF)


imgOut = imgIn( (1: floor(size(imgIn,1)/binF)*binF),  (1: floor(size(imgIn,2)/binF)*binF) );

[m,n]=size(imgOut); 
imgOut = sum( reshape(imgOut,binF,[]) ,1 );
imgOut = reshape(imgOut,m/binF,[]).'; 
imgOut = sum( reshape(imgOut,binF,[]) ,1);
imgOut = reshape(imgOut,n/binF,[]).'; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,fval,exitflag,output] = fminsearchbnd(fun,x0,LB,UB,options,varargin)

xsize = size(x0);
x0 = x0(:);
n=length(x0);

if (nargin<3) || isempty(LB)
  LB = repmat(-inf,n,1);
else
  LB = LB(:);
end
if (nargin<4) || isempty(UB)
  UB = repmat(inf,n,1);
else
  UB = UB(:);
end

if (n~=length(LB)) || (n~=length(UB))
  error 'x0 is incompatible in size with either LB or UB.'
end

if (nargin<5) || isempty(options)
  options = optimset('fminsearch');
end

params.args = varargin;
params.LB = LB;
params.UB = UB;
params.fun = fun;
params.n = n;
params.xsize = xsize;
params.OutputFcn = [];
params.BoundClass = zeros(n,1);
for i=1:n
  k = isfinite(LB(i)) + 2*isfinite(UB(i));
  params.BoundClass(i) = k;
  if (k==3) && (LB(i)==UB(i))
    params.BoundClass(i) = 4;
  end
end

x0u = x0;
k=1;
for i = 1:n
  switch params.BoundClass(i)
    case 1
       if x0(i)<=LB(i)
        x0u(k) = 0;
      else
        x0u(k) = sqrt(x0(i) - LB(i));
      end
   
      k=k+1;
    case 2

      if x0(i)>=UB(i)
          x0u(k) = 0;
      else
        x0u(k) = sqrt(UB(i) - x0(i));
      end
      
       k=k+1;
    case 3
  
      if x0(i)<=LB(i)
       
        x0u(k) = -pi/2;
      elseif x0(i)>=UB(i)
     
        x0u(k) = pi/2;
      else
        x0u(k) = 2*(x0(i) - LB(i))/(UB(i)-LB(i)) - 1;
        x0u(k) = 2*pi+asin(max(-1,min(1,x0u(k))));
      end
      
      k=k+1;
    case 0
      x0u(k) = x0(i);
      
      k=k+1;
    case 4
   
  end
  
end

if k<=n
  x0u(k:n) = [];
end

if isempty(x0u)
  x = xtransform(x0u,params);
  x = reshape(x,xsize);
  fval = feval(params.fun,x,params.args{:});
  exitflag = 0;
  
  output.iterations = 0;
  output.funcCount = 1;
  output.algorithm = 'fminsearch';
  output.message = 'All variables were held fixed by the applied bounds';
  
  return
end

if ~isempty(options.OutputFcn)
  params.OutputFcn = options.OutputFcn;
  options.OutputFcn = @outfun_wrapper;
end

[xu,fval,exitflag,output] = fminsearch(@intrafun,x0u,options,params);
x = xtransform(xu,params);
x = reshape(x,xsize);

  function stop = outfun_wrapper(x,varargin);
    xtrans = xtransform(x,params);
    stop = params.OutputFcn(xtrans,varargin{1:(end-1)});
  end

end 

% ======================================
% ========= begin subfunctions =========
% ======================================
function fval = intrafun(x,params)
xtrans = xtransform(x,params);

fval = feval(params.fun,reshape(xtrans,params.xsize),params.args{:});

end

function xtrans = xtransform(x,params)

xtrans = zeros(params.xsize);

k=1;
for i = 1:params.n
  switch params.BoundClass(i)
    case 1
      xtrans(i) = params.LB(i) + x(k).^2;
      
      k=k+1;
    case 2
      xtrans(i) = params.UB(i) - x(k).^2;
      
      k=k+1;
    case 3
      xtrans(i) = (sin(x(k))+1)/2;
      xtrans(i) = xtrans(i)*(params.UB(i) - params.LB(i)) + params.LB(i);
      xtrans(i) = max(params.LB(i),min(params.UB(i),xtrans(i)));
      
      k=k+1;
    case 4
      xtrans(i) = params.LB(i);
    case 0
      xtrans(i) = x(k);
      
      k=k+1;
  end
end

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mutualInformation = calcMI_v3(im1, im2, qN)


im1uint = int8(qN*im1/max(abs(im1(:))));
im2uint = int8(qN*im2/max(abs(im2(:))));


[~,~,indrow] = unique(im1uint(:)); 
[~,~,indcol] = unique(im2uint(:)); 

jointHistogram = accumarray([indrow indcol], 1);
jointProb = jointHistogram / numel(indrow);
indNoZero = jointHistogram ~= 0;
jointProb1DNoZero = jointProb(indNoZero);
jointEntropy = -sum(jointProb1DNoZero.*log2(jointProb1DNoZero));

histogramImage1 = sum(jointHistogram, 1);
histogramImage2 = sum(jointHistogram, 2);

indNoZero = histogramImage1 ~= 0;
prob1NoZero = histogramImage1(indNoZero);
prob1NoZero = prob1NoZero / sum(prob1NoZero);
entropy1 = -sum(prob1NoZero.*log2(prob1NoZero));

indNoZero = histogramImage2 ~= 0;
prob2NoZero = histogramImage2(indNoZero);
prob2NoZero = prob2NoZero / sum(prob2NoZero);
entropy2 = -sum(prob2NoZero.*log2(prob2NoZero));


mutualInformation = entropy1 + entropy2 - jointEntropy;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%