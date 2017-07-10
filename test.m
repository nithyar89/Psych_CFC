%%Nicola's code


clear
clc

maxNumCompThreads(1);
if ispc
db='C:\Users\nic\Dropbox\'; %LAPTOP
% basedir='';
% rawdir=[basedir 'raw\'];
% matdir=[basedir 'MAT\'];
% scriptdir=[db 'MMN\'];
addpath([db 'dpTOOLS\'])
loadenv;
else 
if exist('/home/nicola/Desktop/Dropbox/','dir'); %LOCAL
    db='/home/nicola/Desktop/Dropbox/';
%     basedir='/home/nicola/Desktop/.../';
elseif exist('/home/polizzotton/Dropbox/','dir') %MEGATRON
    db='/home/polizzotton/Dropbox/';   
%     basedir='/data/data1/users/polizzotton/STCANN_RAW_CONT_OCT2013/';   
    basedir='/data/data1/users/polizzotton/STCANN/'; %June 2104
end
scriptdir=[db 'STCANN/'];
addpath([db 'dpTOOLS/'])
loadenv;
rawdir=[basedir 'Raw_cont/'];
matdir=[basedir 'MAT/'];    
end
addpath(basedir)
addpath(scriptdir)
% 
% runfile_preica=[scriptdir 'RUNpreica.mat'];
% runfile_ica=[scriptdir 'RUNica.mat'];

load chanlocsCNMD
load AREAcnmd

%AVAILABLE FILE NAMES AND NOTES:
% load([scriptdir 'STCANN_SBJINFOwithGR'])
load([scriptdir 'BEHAVBRIEF'])
BEHAVID=ID;
behavid=unique(BEHAVID)';
match=MATCH; acc=ACC; rt=RT; lohi=LOHI;
% [f,id]=occ(ID);
% id=id(f>200)';


FN=dir([rawdir 'S*.raw']); FN={FN.name}';
ID=cell(size(FN));
for n=1:size(FN,1)
   fn=FN{n};   
   ID{n}=fn(strfind(fn,'NN')+2:find(fn=='_',1)-1);
end



gopreica=1;
goica=1;
gopostica=1;
chpostica=1;
gotf=1;

erps=[];


KK=cell(1);
FD=cell(1);
for n=perm(size(ID,1))
try
display([num2str(n), '____________' ID{n} '_________________']) 

fn0=[rawdir FN{n}];
fn1=[matdir 'r' ID{n} '.mat'];
fn1ica=[matdir 'r' ID{n} '_ica.mat']; 

fn2old=[matdir 'fr' ID{n} '.mat']; 
fn2=[matdir 'frnew' ID{n} '.mat']; 
checkfile=[matdir 'checkfile' ID{n} '.mat']; 

% fnTF=[matdir 'tf_' ID{n} '.mat']; %non interpolated

fnTFold=[matdir 'tfint_' ID{n} '.mat']; %interpolated
fnTF_Aold=[matdir 'tfAint_' ID{n} '.mat']; 
fnTF_Bold=[matdir 'tfBint_' ID{n} '.mat']; 
fnTF_ABold=[matdir 'tfABint_' ID{n} '.mat']; 

fnTF=[matdir 'tfnew_' ID{n} '.mat']; %focal comp removal and iterpolated
fnTF_A=[matdir 'tfAnew_' ID{n} '.mat']; 
fnTF_B=[matdir 'tfBnew_' ID{n} '.mat']; 
fnTF_AB=[matdir 'tfABnew_' ID{n} '.mat']; 



if gopreica
if ~exist(fn1,'file') && exist(fn0,'file') 
%CONVERT CONTINUOUS DATA     
[X,sr,~,eegstr]=loadraw2(fn0);

X=fqfiltering(X,sr);    


[epochs,epoch_lookup] =events2epoch('cue1',eegstr,[-.7 3.2],sr);
lo =events2epoch('C1lo',eegstr,[-.7 3.2],sr);
hi =events2epoch('C1hi',eegstr,[-.7 3.2],sr);


cond=zeros(size(epochs(:,1)));
cond(ismember(epochs(:,2),lo(:,2)))=1;
cond(ismember(epochs(:,2),hi(:,2)))=2;

ntr=size(epochs,1);
ntp=1+unique(epochs(:,3)-epochs(:,2));
[v_raw2tr,v_tr2u,v_u2tr]=epochs2vect(epochs); %get reshaping indices
X=reshape(X(:,v_raw2tr),[size(X,1)  ntp ntr]); %3D trialwise format   




ACC=acc(BEHAVID==str2double(ID(n)));
RT=rt(BEHAVID==str2double(ID(n)));
MATCH=match(BEHAVID==str2double(ID(n)));
COND=lohi(BEHAVID==str2double(ID(n)));


%PRE-ICA
MASK=preicaPLUS(X,sr,chanlocs);
X=X(:,v_tr2u);
%Let's store the most compact form (unique copy of trials time points)
display ' '
fprintf(['Saving ' fn1 '...']) 
save(fn1,'X','eegstr','cond','epochs','MASK','ntp','COND','ACC','RT','MATCH')
fprintf(' done.\n')
end
end



if goica
if  exist(fn1,'file') && ~exist(fn1ica,'file')
    display([num2str(n), '____________' ID{n} '_________________']) 
    load(fn1,'X','MASK','epochs')

    
    MASK=MASK(:,:,4);
    btr=all(isnan(MASK));
    bch=all(isnan(MASK),2);
            
    [~,~,~,u2TR,TR2u]=epochs2vect(epochs,[],~btr);
    X=ref(X);
    X=X(:,u2TR);  
    %I need a u2D format    
    X=X(:,TR2u);      
    EYE=findocularsignals(X,sr,chanlocs,find(~bch));      
    X(find(bch),:,:)=[];  
     
%     ieye=[1];
%     [W,L,w,s,time]=decompose(X,EYE(ieye,:));
    
    [W,L,w,s,time]=decompose(X);
   
    fprintf(['Saving ' fn1ica '...']) 
    save(fn1ica,'W','L','w','s','time','btr','bch','EYE')    
    fprintf(' done.\n')    
end
end


if gopostica
if  ~exist(fn2old,'file') && exist(fn1,'file') && exist(fn1ica,'file') 
display([num2str(n), '____________' ID{n} '_____postica________']) 
clear W X   
load(fn1,'X','epochs','MASK','cond')
load(fn1ica,'W','btr','bch','EYE')
MASK=MASK(:,:,4);
[P,i,btr,bch,kTR,salvage]=postica(X,W,MASK,sr,epochs,chanlocs,EYE);
% [P,btr,bch]=justocular(X,MASK,sr,epochs,chanlocs);

display([num2str(n), '____________' ID{n} '___finaldatacleanup_____']) 
[X,btr2,bch,MASK]=finaldatacleanup(P,sr,chanlocs,bch,kTR);
btr(~btr)=btr2;
TR=~btr;
save(fn2old,'X','TR','cond','epochs','MASK','bch','i','salvage')
if exist(fn2old,'file')
    display([fn2old ' was copied!'])
end
end
end


%% NEW
if gopostica
if  ~exist(fn2,'file') && exist(fn1,'file') && exist(fn1ica,'file') && exist(fn2old,'file') 
display([num2str(n), '____________' ID{n} '_____NEWpostica________']) 

load(fn1,'X','epochs','MASK','cond')
load(fn2old,'i','salvage')    
load(fn1ica,'W')

MASK=MASK(:,:,4);  bch=all(isnan(MASK),2);

display([num2str(n), ' FOCAL COMP REMOVAL _ _ _ _']) 
P=reproject(X,W,MASK,sr,epochs,chanlocs,i,1);
display([num2str(n), ' CLEANUP _ _ _ _'])    
[X,btr,bch,MASK]=finaldatacleanup(P,sr,chanlocs,bch);
TR=~btr;

save(fn2,'X','TR','cond','epochs','MASK','bch','i','salvage')
if exist(fn2,'file')
    display([fn2 ' was copied!'])
end
end
end




if  exist(fn2,'file') && ~exist(checkfile,'file')    
    fprintf([num2str(n) ':getting diagnostics...'])     
    load(fn2,'X','TR','cond','MASK')      
    cond=cond(TR);
    
    btr=mean(all(isnan(MASK)));
    
    bch=mean(squeeze(any(isnan(X),2)),2)>.75;
    X(bch,:,:)=NaN;
    
    tr=[sum(cond==1) sum(cond==2)];
    b_tr_ch=[btr mean(bch)];
    
    ERP=erp(X,cond,[],[],[],[1 2]);

    save(checkfile,'ERP','tr','b_tr_ch')
    fprintf(' done.\n')        
end
% 
%     load(checkfile,'ERP')      
%     erps=cat(4,erps,ERP);





if gotf
if  exist(fn2old,'file') && ~exist(fnTFold,'file')    
    load(fn2old,'X','TR','cond')      
    cond=cond(TR);
   
    MASK=squeeze(any(isnan(X),2));
    
% %     %no interpolation 
% %     tr1=sum(~MASK(:,cond==1),2);                 
% %     tr2=sum(~MASK(:,cond==2),2);                 

    %% interpolation
    tr1=ones(size(X,1),1)*sum(cond==1);
    tr2=ones(size(X,1),1)*sum(cond==2);
    X=filleeg(X,chanlocs,'pipe');
    

    display([num2str(n), '_BASIC TIME FREQUENCY ANALYSIS_________'])      
    display('cond==1...')        
    [IND_A,EV_A]=timefq(X(:,:,cond==1),1:70);
    display('cond==2...')        
    [IND_B,EV_B]=timefq(X(:,:,cond==2),1:70);
    save(fnTFold,'IND_*','EV_*','tr1','tr2')
    fprintf(' done.\n')        
end     

if  exist(fn2,'file') && ~exist(fnTF,'file')    
    load(fn2,'X','TR','cond')      
    cond=cond(TR);
   
    MASK=squeeze(any(isnan(X),2));
    
% %     %no interpolation 
% %     tr1=sum(~MASK(:,cond==1),2);                 
% %     tr2=sum(~MASK(:,cond==2),2);                 

    %% interpolation
    tr1=ones(size(X,1),1)*sum(cond==1);
    tr2=ones(size(X,1),1)*sum(cond==2);
    X=filleeg(X,chanlocs,'pipe');
    

    display([num2str(n), '_BASIC TIME FREQUENCY ANALYSIS_________'])      
    display('cond==1...')        
    [IND_A,EV_A]=timefq(X(:,:,cond==1),1:70);
    display('cond==2...')        
    [IND_B,EV_B]=timefq(X(:,:,cond==2),1:70);
    save(fnTF,'IND_*','EV_*','tr1','tr2')
    fprintf(' done.\n')        
end
end



if  exist(fnTFold,'file') && ~exist(fnTF_ABold,'file')  
    display([num2str(n), '_COLLAPSING TIME FREQUENCY_________'])   
    
    load(fnTFold,'IND_*','tr1','tr2')
    ind_A=IND_A; ind_B=IND_B;
    ind_A(:,:,tr1<50)=NaN; %at least 50 trials per electrode?
    ind_B(:,:,tr2<50)=NaN;
    
    fprintf('2AREA... ')   
    IND_A=el2a(ind_A,AREA);
    IND_B=el2a(ind_B,AREA);
    save(fnTF_Aold,'IND_*')
    fprintf(' done.\n')    
    
    
    fprintf('2BAND... ')  
    IND_A=f2b(ind_A,1:size(ind_A,1));
    IND_B=f2b(ind_B,1:size(ind_B,1)); 
    save(fnTF_Bold,'IND_*')
    fprintf(' done.\n')    
    
 
    fprintf('2AB... ')   
    IND_A=el2a(f2b(ind_A,1:size(ind_A,1)),AREA);
    IND_B=el2a(f2b(ind_B,1:size(ind_B,1)),AREA);
    save(fnTF_ABold,'IND_*')
    fprintf(' done.\n')        
end
if  exist(fnTF,'file') && ~exist(fnTF_AB,'file')  
    display([num2str(n), '_COLLAPSING TIME FREQUENCY_________'])   
    
    load(fnTF,'IND_*','tr1','tr2')
    ind_A=IND_A; ind_B=IND_B;
    ind_A(:,:,tr1<50)=NaN; %at least 50 trials per electrode?
    ind_B(:,:,tr2<50)=NaN;
    
    fprintf('2AREA... ')   
    IND_A=el2a(ind_A,AREA);
    IND_B=el2a(ind_B,AREA);
    save(fnTF_A,'IND_*')
    fprintf(' done.\n')    
    
    
    fprintf('2BAND... ')  
    IND_A=f2b(ind_A,1:size(ind_A,1));
    IND_B=f2b(ind_B,1:size(ind_B,1)); 
    save(fnTF_B,'IND_*')
    fprintf(' done.\n')    
    
 
    fprintf('2AB... ')   
    IND_A=el2a(f2b(ind_A,1:size(ind_A,1)),AREA);
    IND_B=el2a(f2b(ind_B,1:size(ind_B,1)),AREA);
    save(fnTF_AB,'IND_*')
    fprintf(' done.\n')        
end


catch err
display('******ERROR*******')    
display(err.message)
end
end






break




%% STATS
clear

rmpath(genpath('/raid5/rcho/TOOLS/NIC/MATLAB/fieldtrip-20150503/'));%removing overlapping functions
rmpath(genpath('/raid5/rcho/TOOLS/NIC/MATLAB/eeglab13_4_4b/functions/octavefunc/'));%removing overlapping functions
if ispc
db='C:\Users\nic\Dropbox\'; %LAPTOP
% basedir='C:\Users\nic\Desktop\\';
addpath([db 'dpTOOLS\'])
loadenv;
else
if exist('/raid5/rcho/STCANN/STCANN/','dir'); %GROND
db='/raid5/rcho/STCANN/STCANN/';
basedir='/raid5/rcho/STCANN/swm/adults/';
end
addpath([db])% addpath([db 'dpTOOLS/'])
loadenv;
end
scriptdir=[db];% scriptdir=[db 'STCANN/'];
addpath(scriptdir)
load chanlocsCNMD
load AREAcnmd2
sr=250;
rawdir=[basedir '/Raw_cont/'];
matdir=[basedir '/MAT/'];
addpath(basedir)
addpath(scriptdir)
addpath(matdir)

stataname='STATANEWEST';
KEY{1}='tfABint_';
KEY{2}='tfAint_';
KEY{3}='tfBint_';
KEY{4}='tfint_';

t=(0:976-1)*4; t=t-700;
bl=t>-500 & t <-150;

it=t>-100&t<2100;


outlierremoval=0;
%nk=2
for nk=1:3
key=KEY{nk};
fn=[matdir key '*.mat'];
file=dir(fn); file={file.name}';
clear id
file=sort(file);
for n=1:size(file,1)
    a=file{n};
    id(n,1)=str2double(a((find(a=='_',1,'last')+1):end-4));
end

% % badsbj=ismember(id,[671 603 766]);
% % id(badsbj)=[];

load ELN






load([matdir file{1}])
dims=size(IND_A);

A=nan([dims size(file,1)]);
B=nan([dims size(file,1)]);
if outlierremoval
fn2=[matdir stataname '_corrected_' key '.mat'];   
else
fn2=['/raid5/rcho/STCANN/swm/adults/' stataname '_' key '.mat'];
%fn2=[matdir stataname '_' key '.mat'];       
end

display(key)
for n=1:size(file,1)
display([key ': gathering data for ' file{n}])            
load([matdir file{n}],'IND_A','IND_B')
A(:,:,:,n)=IND_A;
B(:,:,:,n)=IND_B;
end

bads=[6    12    24    31    56];
A(:,:,:,bads)=[];
B(:,:,:,bads)=[];
id(bads)=[];

GR=makeas(IDELN,id,ELN);


if outlierremoval
    
A=baseline(A,bl,2); 
B=baseline(B,bl,2); 


% A=A(:,it,:,:);
% B=B(:,it,:,:);
%%optional%%%%%%%%%%%%%%%%%%%%%%%%%
a=squeeze(all(squeeze(all(isnan(A),2)),1));
b=squeeze(all(squeeze(all(isnan(B),2)),1));
badsubject=mean(a)>.3 | mean(b)>.3;
A=A(:,:,:,~badsubject);
B=B(:,:,:,~badsubject);
file=file(~badsubject);
GR=GR(~badsubject);
TT2=[];
for nl=3:size(A,1)
display(['correcting data for layer ' num2str(nl)])            
    
LL=zeros([size(A,3),size(A,2),2,size(A,4)]);
for n=1:size(A,4)
LL(:,:,1,n)=squeeze(A(nl,:,:,n))';
LL(:,:,2,n)=squeeze(B(nl,:,:,n))';
end   
LL2=LL-repmat(nanmean(LL(:,:),2),[1 size(LL,2) size(LL,3) size(LL,4) ]);
LL2=abs(LL2./repmat(nanstd(LL2(:,:),1,2),[1 size(LL,2) size(LL,3) size(LL,4) ]));
LL2=repmat(any(LL2>10,2),[1 size(LL2,2) 1 1]);
LL(LL2)=NaN;    
  
[el tp cond sbj]=size(LL);
LL2=zeros([el cond*sbj*tp]);
k=0;
for s=1:sbj
for c=1:cond    
k=k+1;   
x=LL(:,:,c,s);
LL2(:,tp*(k-1)+(1:tp))=x;
end
end 

gel=~any(isnan(LL2),2);
if any(~gel)
fprintf('Filling missing values with ppca ...')
rk=[rank(LL2(gel,:)) round(size(LL2,1)/2)];
[~, ~, ~, ~,LL2adj] = ppca_mv(LL2',max(rk),0);
LL2adj=LL2adj';
LL2(isnan(LL2))=LL2adj(isnan(LL2));
clear LL2adj
fprintf('done.\n')
end



[~,~,E,T2]=princomp(LL2');
% [L,PC,E,T2]=princomp(LL2');
% PC=PC';
E=E/sum(E);
c2retain=find(cumsum(E)>.9,1,'first')-1;

 
k=0;
t2=LL(1,:,:,:,:);
for s=1:sbj
for c=1:cond    
k=k+1;   
LL(:,:,c,s)=LL2(:,tp*(k-1)+(1:tp));
t2(:,:,c,s)=T2(tp*(k-1)+(1:tp));
% PPCC(:,:,c,s)=PC(:,tp*(k-1)+(1:tp));
end
end
TT2(nl,:)=max(meandim(t2,2));

% [37 24 84] i.e id=[671 603 766];

for n=1:size(A,4)
A(nl,:,:,n)=reshape(LL(:,:,1,n)',[1 size(LL(:,:,1,n)')] );
B(nl,:,:,n)=reshape(LL(:,:,2,n)',[1 size(LL(:,:,1,n)')] );
end  
end


TT2=TT2';
TT2=TT2(:,3:end-1);
badsubject=any(zscore(TT2)>3,2);
A=A(:,:,:,~badsubject);
B=B(:,:,:,~badsubject);
GR=GR(~badsubject);
A=nanwfilter(A,20,2);
B=nanwfilter(B,20,2);
[T,M,P,TG,MG,PG,bgT,bgP]=tfstats(A,B,GR);
t=(0:676-1)*4; t=t-600;
t=t(it);
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    keyboard
nanA=isnan(A); nanB=isnan(B); %this section was added around the wfilter call b/c the NANs were preventing filtfilt from running properly
A(nanA)=0; B(nanB)=0;
A=nanwfilter(A,20,2);  %this function doesn't accept NANs so it was removed
B=nanwfilter(B,20,2);
A(nanA)=nan; B(nanB)=nan;


if nk==4
    A=baseline(A,bl);
    B=baseline(B,bl);
A=A(:,5:5:end,:,:);
B=B(:,5:5:end,:,:);
[T,M,P,TG,MG,PG,bgT,bgP]=tfstats(A,B,GR);    
else
[T,M,P,TG,MG,PG,bgT,bgP]=tfstats(A,B,GR,bl);
end
end

save(fn2,'T','M','P','TG','MG','PG','bgT','bgP','GR','t')
end


break


imagx(-T,3,t,[],[0 300 500 800 1000 1300 2800])


GRlbl=[{'EARLY'} {'LATE'} {'NO'}];
GRlbl=[{'EARLY-LATE'} {'EARLY-NO'} {'LATE-NO'}];

for n=1:8
for g=1:3    
LBL{n,g}=[GRlbl{g} '-' AREA(n).name]
end
end





it=t>-400&t<2000;
BB=zeros([129,sum(it),2,size(A,4)]);
for n=1:size(A,4)
BB(:,:,1,n)=squeeze(A(6,it,:,n))';
BB(:,:,2,n)=squeeze(B(6,it,:,n))';
end


[el tp cond sbj]=size(BB);
BB2=zeros([el cond*sbj*tp]);
k=0;
for s=1:sbj
for c=1:cond    
k=k+1;   
x=BB(:,:,c,s)-repmat(mean(BB(:,1:100,c,s),2),[1 size(BB,2)]);
BB2(:,tp*(k-1)+(1:tp))=x;
end
end 

gel=~any(isnan(BB2),2);
if any(~gel)
fprintf('Filling missing values with ppca ...')
rk=[rank(BB2(gel,:)) round(size(BB2,1)/2)];
[~, ~, ~, ~,BB2adj] = ppca_mv(BB2',max(rk),0);
BB2adj=BB2adj';
BB2(isnan(BB2))=BB2adj(isnan(BB2));
fprintf(' done.\n')
end


[iW,C,E,T2]=princomp(BB2');
C=C';
E=E/sum(E);
c2retain=find(cumsum(E)>.9,1,'first')-1;


 
 
k=0;
CC=BB(1:size(C,1),:,:,:);
t2=BB(1,:,:,:,:);
for s=1:sbj
for c=1:cond    
k=k+1;   
CC(:,:,c,s)=C(:,tp*(k-1)+(1:tp));
t2(:,:,c,s)=T2(tp*(k-1)+(1:tp));
end
end
T2=t2;
t2=meandim(T2,2);
T2=squeeze(T2);
clear C
 
 
outliers=any(t2<1) | any(transf(t2)>3)
d=diff(t2)./sum(t2);
d=zscore(d(~outliers));

 
 
 ACC=[];
 RT=[];
 ID2=[];
 COND=[];
 for n=1:size(S,2)
   ID2=[ID2 S(n).ID*ones(size(S(n).ACC))];  
   ACC=[ACC S(n).ACC]; 
   COND=[COND S(n).COND(1:numel(S(n).ACC))];    
   RT=[RT S(n).RT]; 
 end
 
 CAT=[COND' ACC' RT'];
 
 
 
 
 
 
 
 

 
 
 
 
 
