function [X,lambda,it] = basic(A,p,X_0,itmax,tol)
%function [X,lambda,it] = basic(A,p,X_0,itmax,tol)
%Esegue il Metodo delle Iterazioni Ortogonali nella versione standard, determinando le p autocoppie dominanti della matrice A a partire dalla matrice con colonne ortonormali X_0. 
%IN INPUT
%   - p: numero di autocoppie dominanti da calcolare
%   - X_0: matrice con colonne ortonormali di partenza
%   - itmax: numero massimo di iterazioni svolte
%   - tol: tolleranza per condizioni di terminazione
%IN OUTPUT
%   - X: approssimazioni degli autovettori 
%   - lambda: approssimazioni degli autovalori
%   - it: numero di iterazioni effettivamente svolte

    X = X_0;
    n = size(X,1);
    lambda = zeros(p,1);
    for it = 1:itmax
        lambda_prec = lambda;            %precedente, necessario per la condizione di terminazione
        X = A*X;
        [X,R] = qr(X,0);
        lambda = diag(R);
        if (norm(lambda_prec-lambda)/norm(lambda) < tol)
            break;
        end
    end
end