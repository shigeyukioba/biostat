%% Loading example data
% ���K�p�̃f�[�^��MATLAB�̒��ɓǂݍ���

X_org = load('NBLexpression.dat');

z = X_org(1,:);
z
X = X_org(2:end,:);
size(X)

%% Exercise 1. HISTGRAM
% Here, we draw a histgram of expression levels of a gene with specifying mean and SD.
% ��̈�`�q�̔����ʂ̃q�X�g�O�����i�x�����z�j����������B���ςƕW���΍������킹�ĕ\������B

figure
i1 = 1; gname1 = sprintf('gene%d',i1);
hist( X(i1, :) )
xlabel( gname1 )
% Add mean and SD.
m = mean( X(i1,:) )
sd = std( X(i1,:) )
hold on
plot( m+[0,0], [0,25], 'r-', 'LineWidth', 2 )
plot( m+sd+[0,0], [0,25], 'r--', 'LineWidth', 2 )
plot( m-sd+[0,0], [0,25], 'r--', 'LineWidth', 2 )

%% You can draw multiple histgram by using subplot
% subplot ��p���邱�Ƃŕ����̓x�����z��`���Ă݂܂�

figure
i0 = [1, 3, 5, 10, 12, 15];
for i = 1:length(i0)
   subplot(2,3,i)
   hist( X(i0(i),:) )
   gname = sprintf('gene%d',i0(i));
   xlabel( gname )
end

%% Exercise 2. SCATTER PLOT
% �Q�̈�`�q�̔����ʎU�z�}�iscatter plot�j��`���܂��傤

i1 = 2; i2 = 3;
figure
plot( X(i1,:),X(i2,:),'.')
xlabel( sprintf('gene%d',i1) )
ylabel( sprintf('gene%d',i2) )


%% Specify the MEAN
% ���ϒl�����킹�ĕ\�����܂�

m1 = mean( X(i1,:) );
m2 = mean( X(i2,:) );
hold on
plot( m1, m2, 'ro' )


%% Add a unit circle
% �P�ʉ~��`���Ă݂܂��傤

t = 2*pi*([0:360]/360);
gx = sin(t);
gy = cos(t);
plot( gx+m1, gy+m2, 'g-' )

%% Calculate standard deviation
% �W���΍����v�Z���܂�
sd1 = std( X(i1, :) )
sd2 = std( X(i2, :) )

%plot( sd1*gx+m1, sd2*gy+m2, 'b-' )

%% Calculate variance-covariance-matrix ���U�����U�s�� C ���v�Z���܂�
C = cov( X([i1, i2], :)' )

%% Calculate eigenvalue decomposition of C. �s�� C �̌ŗL�l�������v�Z���܂� 
[V,D] = eig( C )

% Be sure that  C = V * D * V'  holds. �ŗL�l�����̌��ʂ��m���߂܂�
sum( sum( abs( V*D*V' - C ) ) )

%% Draw a covariance oval
% �����U�ȉ~��`���܂�

ex = V(1,1)*sqrt(D(1,1))*gx + V(1,2)*sqrt(D(2,2))*gy;
ey = V(2,1)*sqrt(D(1,1))*gx + V(2,2)*sqrt(D(2,2))*gy;
plot( ex+m1, ey+m2, 'r-', 'LineWidth', 2 )

%% Exercise for Report: 3 multiple scatter plots
% Draw six scatter plots with covariance ovals for six arbitrary pairs of
% genes.
% ��`�q i1 �ƈ�`�q i2 �̃y�A�ɑΉ�����U�z�}���ȉ��ŕ\�����邱�Ƃ��ł��܂��B
% ���̂����A�\��ǁiz=0�j�Ɨ\�㈫�iz=1�j�̗��ƐԂŐF�������܂��B
% ������Q�l�ɂ��ĔC�ӂ�6�y�A�ɑΉ�����U�z�}�������U�ȉ~�����ŕ\�����Ă݂܂��傤�B
% 6�̈�`�q�y�A�̒��ɗ\��\���ɖ𗧂������ȃy�A�͂���܂������H

figure
plot( X(i1,z==0), X(i2,z==0), 'b.', ...
    X(i1,z==1),X(i2,z==1), 'rx')
xlabel( sprintf('gene%d',i1) )
ylabel( sprintf('gene%d',i2) )


%% Exercise 4. PCA
% �听�����͂��s�����ƂŁA���听���Ƒ��听�������߁A�U�z�}��`���܂��B
% ���̂����A�\��ǁiz=0�j�Ɨ\�㈫�iz=1�j�̗��F�������܂��B

geneset = [1:500]; % an arbitrary geneset
Xo = X( geneset, : );
m = mean( Xo, 2 ); % m(i) is mean expr. level of gene  i
Y = Xo - repmat( m, 1, size(Xo,2) );
[U,S,V] = svds( Y, 30 ); % singular value decomposition
size(V)

figure
plot( V(z==0,1), V(z==0,2), 'b.', ...
    V(z==1,1), V(z==1,2), 'ro')
xlabel('PC1')
ylabel('PC2')


