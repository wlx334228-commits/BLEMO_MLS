classdef CGANv9
    %该版本的最终判别器是通过每若干代训练所得的弱判别器集成所得，集成的权重由弱判别器对应的损失函数确定。            
    %   此处显示详细说明
    
    properties
        l_upper;
        l_lower;
        
        latent_dim;
        numnodes_G;
        numnodes_D;
        
        paramsGen;
        paramsDis;
        stGen;
        stDis;
        avgG;
        avgGS;
        
        basis_paramsDis;
        D_loss;
        m;
        w;
    end
    
    methods
        function obj = CGANv9(UX,LX,Global)
            %UNTITLED2 构造此类的实例
            %   此处显示详细说明
            obj.l_upper = Global.lower_domain(2,:);
            obj.l_lower = Global.lower_domain(1,:);
            
            %% Settings
            maxepochs = Global.D(2)*20;
            
            obj.latent_dim = Global.D(2);
            obj.numnodes_G = Global.D(2)*3; %size(LX,2)+
            obj.numnodes_D = Global.D(2)*3+2;%size(UX,2)+size(LX,2)+
            %% Initialization
            %% Generator
            paramsGen.EMW1 = dlarray(...
                initializeGaussian([obj.numnodes_G,...
                obj.latent_dim],1));
            paramsGen.EMb1 = dlarray(zeros(obj.numnodes_G,1,'double'));
            
            paramsGen.FCW1 = dlarray(...
                 initializeGaussian([obj.numnodes_G,obj.numnodes_G+Global.D(1)],1));            
            paramsGen.FCb1 = dlarray(zeros(obj.numnodes_G,1,'double'));            
            paramsGen.BNo1 = dlarray(zeros(obj.numnodes_G,1,'double'));
            paramsGen.BNs1 = dlarray(ones(obj.numnodes_G,1,'double'));
            paramsGen.FCW2 = dlarray(initializeGaussian([obj.numnodes_G,obj.numnodes_G+Global.D(1)],1));
            paramsGen.FCb2 = dlarray(zeros(obj.numnodes_G,1,'double'));
            paramsGen.BNo2 = dlarray(zeros(obj.numnodes_G,1,'double'));
            paramsGen.BNs2 = dlarray(ones(obj.numnodes_G,1,'double'));
            
            paramsGen.FCW3 = dlarray(initializeGaussian([Global.D(2),obj.numnodes_G+Global.D(1)],1));
            paramsGen.FCb3 = dlarray(zeros(Global.D(2),1,'double'));            
            
            stGen.BN1 = []; stGen.BN2 = []; 
            
            obj.paramsGen = paramsGen;   obj.stGen = stGen;
            %% Discriminator
            paramsDis.FCW1 = dlarray(initializeGaussian([obj.numnodes_D,...
                Global.D(2)],1));          
            paramsDis.FCb1 = dlarray(zeros(obj.numnodes_D,1,'double'));                      
            paramsDis.BNo1 = dlarray(zeros(obj.numnodes_D,1,'double'));
            paramsDis.BNs1 = dlarray(ones(obj.numnodes_D,1,'double'));
            paramsDis.FCW2 = dlarray(initializeGaussian([obj.numnodes_D,obj.numnodes_D]));
            paramsDis.FCb2 = dlarray(zeros(obj.numnodes_D,1,'double'));
            paramsDis.BNo2 = dlarray(zeros(obj.numnodes_D,1,'double'));
            paramsDis.BNs2 = dlarray(ones(obj.numnodes_D,1,'double'));
            paramsDis.FCW3 = dlarray(initializeGaussian([1,obj.numnodes_D]));
            paramsDis.FCb3 = dlarray(zeros(1,1,'double'));
            
            paramsDis.EMW1 = dlarray(...
                 initializeUniform([obj.numnodes_D,...
                 Global.D(1)]));
            paramsDis.EMb1 = dlarray(zeros(obj.numnodes_D,1,'double')); 
            stDis.BN1 = []; stDis.BN2 = [];
            
            % average Gradient and average Gradient squared holders
            avgG.Dis = []; avgGS.Dis = []; avgG.Gen = []; avgGS.Gen = [];
            
            obj.paramsDis = paramsDis; obj.stDis = stDis; obj.avgG = avgG; obj.avgGS = avgGS;
            
            obj.basis_paramsDis = [];
            obj.D_loss = [];
            obj.w  = [];
            
            %% Train
            obj = train(obj,UX,LX,maxepochs);
        end
        
        
        function obj = train(obj,UX,LX,maxepochs)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            
            %% Settings
            batch_size = 10;
            lrD = 0.004;
            lrG = 0.005;
            beta1 = 0.5;
            beta2 = 0.999;
            k = 1;
            
                                  
            obj.m = floor(maxepochs/10);
            %% Load data
            [n,K] = size(UX);                 
            trainX = (LX)';
            trainY = (UX)';
                        
%             Xt1 = UX;
           
            %% Train
            
            numIterations = floor(n/batch_size);
            out = false; epoch = 0; global_iter = 0;           
            while ~out
                
%                 if  mod(epoch+1,10)==0
%                     
%                    if length(obj.basis_paramsDis)<obj.m
%                       obj.basis_paramsDis = [obj.basis_paramsDis,obj.paramsDis];
%                       obj.D_loss = [obj.D_loss,d_loss];                     
%                    else
%                        [~,I] = max(obj.D_loss);
%                        if d_loss <=obj.D_loss(I)
%                            obj.basis_paramsDis(I) = obj.paramsDis;
%                            obj.D_loss(I) = d_loss; 
%                        end
%                    end
%                    
%                 end
                
%                  xt = [rand*3+1,rand(1,K-1)*2*K-K];
%                  Xt2 = repmat(xt,50,1);
                 shuffleid = zeros(k+1,n);
                 for i = 1:k+1
                    shuffleid(i,:) = randperm(n);%打乱训练样本顺序
                 end

                for i=1:numIterations
                    global_iter = global_iter+1;                    
                    idx = (i-1)*batch_size+1:i*batch_size;
                
                    
                    for j = 1:k
                       
                        XBatch=gpdl(double(trainX(:,shuffleid(j,idx))),'CB');      %小批次训练样本
                        YBatch=gpdl(double(trainY(:,shuffleid(j,idx))),'CB');      %小批次训练样本的标签
                        
                        
                        noise = gpdl(randn([obj.latent_dim,...     %产生服从正态分布的噪声
                            batch_size]),'CB');
                        
                        [GradDis,obj.stGen,obj.stDis,d_loss] = ...                  %计算误差梯度
                            dlfeval(@modelGradients_D,XBatch,YBatch,noise,...
                            obj.paramsGen,obj.paramsDis,obj.stGen,obj.stDis);                       
                        
                        % Update Discriminator network parameters
                        [obj.paramsDis,obj.avgG.Dis,obj.avgGS.Dis] = ...                 %根据误差梯度更新判别器
                            adamupdate(obj.paramsDis, GradDis, ...
                            obj.avgG.Dis, obj.avgGS.Dis, global_iter, ...
                            lrD, beta1, beta2);
                    end
                    
                    noise = gpdl(randn([obj.latent_dim,...     %产生服从正态分布的噪声
                            batch_size]),'CB');
                    
                    YBatch=gpdl(double(trainY(:,shuffleid(k+1,idx))),'CB');         %小批次训练样本的标签
                    
                    GradGen = ...                  %计算误差梯度
                            dlfeval(@modelGradients_G,YBatch,noise,...
                            obj.paramsGen,obj.paramsDis,obj.stGen,obj.stDis);       
                    
                    % Update Generator network parameters
                    [obj.paramsGen,obj.avgG.Gen,obj.avgGS.Gen] = ...                %根据误差梯度更新生成器
                        adamupdate(obj.paramsGen, GradGen, ...
                        obj.avgG.Gen, obj.avgGS.Gen, global_iter, ...
                        lrG, beta1, beta2);
                end              
                
%-----------------------------Test——————————————————————————————————————————               
                
%                 Yt1 = obj.Generator(Xt1);
%                 Yt2 = obj.Generator(Xt2);
%                 p = obj.Discriminator(Xt1,Yt1);
%                 I = find(p<0);
%                 
%                 subplot(1,2,1)
%                 Draw(LX,'rs')
%                 Draw(UX,'ro');
%                 Draw(Yt1(I,:),'bs');
%                 Yt1(I,:) = [];
%                 Draw(Yt1,'cs');
%                 Draw(Xt1,'bo');
% 
%                 Yp = rand(1000,K)*2*K-K;
%                 p = obj.Discriminator(repmat(xt,1000,1),Yp);
%                 
%                 I = find(p<=0.2);
%                 
%                 subplot(1,2,2)                
%                 Draw(Yp(I,:),'c.');
%                 Yp(I,:) = [];
%                 Draw(Yp,'r.');
%                 Draw(Yt2,'bs');
%                 Draw(Xt2,'ro');
% 
%                 hold off;
%                 pause(0.001);
%---------------------------------------------------------------------------------------                
                
                epoch = epoch+1;
                if epoch == maxepochs
                    out = true;
                end 
                
            end
%             
%             max_Dloss = max(obj.D_loss);
%             obj.w = (max_Dloss - obj.D_loss)/sum(max_Dloss - obj.D_loss);


        end
        function LX = Generator(obj,UX)
            n = size(UX,1);
            
            testY = gpdl(double(UX'),'CB');
            noise = gpdl(randn([obj.latent_dim,...     %产生服从正态分布的噪声
                            n]),'CB');
            testX = generator(noise,testY,obj.paramsGen,obj.stGen);
            
            LX = extractdata(stripdims(testX))';
            
            %修正
            l_uppers = repmat(obj.l_upper,n,1);
            l_lowers = repmat(obj.l_lower,n,1);
            LX = min(max(LX,l_lowers),l_uppers);
        end
        
        function p = Discriminator(obj,UX,LX)
            testY = gpdl(double(UX'),'CB');
            testX = gpdl(double(LX'),'CB');
            params = obj.paramsDis;
            

            %1
            dly = fullyconnect(testX,params.FCW1,params.FCb1);
%            dly = leakyrelu(dly,0.2);
            dly = elu(dly,1);
            %2
            dly1 = fullyconnect(dly,params.FCW2,params.FCb2);
%             dly = leakyrelu(dly,0.2);
            dly1 = elu(dly1,1);
%             %3
            dly1 = fullyconnect(dly1,params.FCW3,params.FCb3);
            % sigmoid
%            dly = sigmoid(dly);
            dly2 = embedding(dly,testY,params);
            
            dly = dly1+dly2;

           p = extractdata(stripdims(dly))';
        end
        
         function p = ensemble_Discriminator(obj,UX,LX)
            testY = gpdl(double(UX'),'CB');
            testX = gpdl(double(LX'),'CB');
            
            P = zeros(1,size(testY,2));
            
            for i = 1:obj.m
                params = obj.basis_paramsDis(i);

                % fully connected
                %1
                dly = fullyconnect(testX,params.FCW1,params.FCb1);
%               dly = leakyrelu(dly,0.2);
                dly = elu(dly,1);
                %2
                dly1 = fullyconnect(dly,params.FCW2,params.FCb2);
%               dly = leakyrelu(dly,0.2);
                dly1 = elu(dly1,1);
%               %3
                dly1 = fullyconnect(dly1,params.FCW3,params.FCb3);
                % sigmoid
%               dly = sigmoid(dly);
                dly2 = embedding(dly,testY,params);
            
                dly = dly1+dly2;
            
                P = P + obj.w(i)*dly;
            end
            p = extractdata(stripdims(P))';
        end
    end
end

%% Helper Functions
%% gpu dl array wrapper
function dlx = gpdl(x,labels)
% dlx = gpuArray(dlarray(x,labels));
dlx = dlarray(x,labels);
end
%% Weight initialization
function parameter = initializeGaussian(parameterSize,sigma)
if nargin < 2
    sigma = 0.05;
end
parameter = randn(parameterSize, 'double') .* sigma;
end
function parameter = initializeUniform(parameterSize,sigma)
if nargin < 2
    sigma = 0.05;
end
parameter = 2*sigma*rand(parameterSize, 'double')-sigma;
end
%% Generator
function [dly,st] = generator(dlx,labels,params,st)

dly = fullyconnect(dlx,params.EMW1,params.EMb1);
% dly = elu(dly,1);
dly = [dly;labels];
%1
dly = fullyconnect(dly,params.FCW1,params.FCb1);   %dly = params.FCW1 * dly + params.FCb1
%dly = leakyrelu(dly,0.2);                          %激活层
dly = elu(dly,1);
if isempty(st.BN1)
    [dly,st.BN1.mu,st.BN1.sig] = batchnorm(dly,...
        params.BNo1,params.BNs1,'MeanDecay',.8);
else
    [dly,st.BN1.mu,st.BN1.sig] = batchnorm(dly,params.BNo1,...
        params.BNs1,st.BN1.mu,st.BN1.sig,...
        'MeanDecay',.8);
end
dly = [dly;labels];
%2
dly = fullyconnect(dly,params.FCW2,params.FCb2);
%dly = leakyrelu(dly,0.2);
dly = elu(dly,1);
if isempty(st.BN2)
    [dly,st.BN2.mu,st.BN2.sig] = batchnorm(dly,...
        params.BNo2,params.BNs2,'MeanDecay',.8);
else
    [dly,st.BN2.mu,st.BN2.sig] = batchnorm(dly,params.BNo2,...
        params.BNs2,st.BN2.mu,st.BN2.sig,...
        'MeanDecay',.8);
end
dly = [dly;labels];
 dly = fullyconnect(dly,params.FCW3,params.FCb3);

end
%% Discriminator
function [dly,st] = discriminator(dlx,labels,params,st)

%1
dly = fullyconnect(dlx,params.FCW1,params.FCb1);
%dly = leakyrelu(dly,0.2);
dly = elu(dly,1);
dly = dropout(dly,.3);
%2
dly1 = fullyconnect(dly,params.FCW2,params.FCb2);
%dly = leakyrelu(dly,0.2);
dly1 = elu(dly1,1);
dly1 = dropout(dly1,.3);
% 3
dly1 = fullyconnect(dly1,params.FCW3,params.FCb3);
% sigmoid
%dly = sigmoid(dly);

dly2 = embedding(dly,labels,params);

dly = dly1 + dly2;
end
%% modelGradients
function [GradDis,stGen,stDis,d_loss]=modelGradients_D(x,y,z,paramsGen,...
    paramsDis,stGen,stDis)

%用于计算判别器误差
[fake_images,stGen] = generator(z,y,paramsGen,stGen);                   %生成器生成一批次假样本
d_output_real = discriminator(x,y,paramsDis,stDis);                     %判别器判别一批次真样本

[d_output_G,stDis] = discriminator(fake_images,y,paramsDis,stDis);   %判别器判别与真样本具有相同标签的假样本

% Loss due to true or not
d_loss =mean(max(0,1-d_output_real)+ max(0,1+d_output_G));

% For each network, calculate the gradients with respect to the loss.
GradDis = dlgradient(d_loss,paramsDis);

d_loss = extractdata(stripdims(d_loss));
end

function GradGen=modelGradients_G(y,z,paramsGen,...
    paramsDis,stGen,stDis)


%用于计算生成器误差
fake_images0 = generator(z,y,paramsGen,stGen);                         %判别器判别另一批次假样本
d_out_fake0 = discriminator(fake_images0,y,paramsDis,stDis);

%g_loss = -mean(log(d_out_fake0+ eps));
g_loss = -mean(d_out_fake0);
% vector = (1-2*round(extractdata(d_out_fake0)));
% SL = ones(size(vector));
% SL(vector<0)=0.6;
% g_loss = -mean(SL.*vector.*log(d_out_fake0+ eps));


% For each network, calculate the gradients with respect to the loss.
% GradGen = dlgradient(g_loss,paramsGen,'RetainData',true);
GradGen = dlgradient(g_loss,paramsGen);
end

%% dropout
function dly = dropout(dlx,p)
if nargin < 2
    p = .3;
end
[n,d] = rat(p);
mask = randi([1,d],size(dlx));
mask(mask<=n)=0;
mask(mask>n)=1;
dly = dlx.*mask;

end
%% embedding
function dly = embedding(dlx,labels,params)

dly = fullyconnect(labels,params.EMW1,params.EMb1);
% dly = leakyrelu(dly,0.2);
dly = elu(dly,1);

% dly = fullyconnect(dly,params.EMW2,params.EMb2);
% % dly = leakyrelu(dly,0.2);
% dly = elu(dly,1);

dly = sum(dly.*dlx);


end

function y = elu(x,alpha)
    y = max(0,x)+ min(0,alpha*(exp(x)-1));
end