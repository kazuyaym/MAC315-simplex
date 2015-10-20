% EP2 - MAC0315 Programação Linear
% Marcos Kazuya Yamazaki 7577622
% Victor de Queiros Mattoso Archela dos Santos 7557202

function [ind v] = simplex(A,b,c,m,n,x)
	ind = 1; % ind vale -1 caso o PL eh ilimitado, ou 0 caso tenha sol otima
	v = x; % caso o PL tenha solucao otima, v eh a sol otima
	
	iterando = 0; % Numero de iteracoes do metodo simplex revisado
	
	printf("Simplex: Fase 2\n\n");

    IB = zeros (0, 1); % indices das variaveis basicas    
    IN = zeros (0, 1); % indices das variaveis nao-basicas
    
    % matriz A somente com as colunas basicas, que no comeco
    % ela eh uma matriz nao inversa, mas com o decorrer das iteracoes 
    % somente a matriz inversa relacionadas a base sera mantida, pois
    % a matriz basica nao eh usada
    B  = zeros (m, 0);
    
    % ----------------------------------------------------------------- PASSO 0:
    %
    % Antes do inicio da iteracao conhecemos uma base formada 
    % pelas colunas da matriz A, das variaveis basicas.
    % Assim montamos a matriz B com essas colunas.
    % E tambem conhecemos a matriz inversa de B = inv(B)

    for i = 1:n
        if ( x(i) > 0 )
		    IB = [IB; i];
		    B  = [ B, A(:,i)];
	    else
	        IN = [IN; i];
	    endif
    endfor
    
    B = inv(B); % Matriz inversa calculada apenas uma vez antes da iteracao
	    
	while(ind == 1)
	    printf("----------- Iterando %2d -----------\n", iterando++);
	    for i = 1:m
	        printf("%d %.5f\n", IB(i), x(IB(i)));
	    endfor
	    printf("\n");
	    printf("Valor da função objetivo: %.5f\n\n", c'*x);

	    % ------------------------------------------------------------- PASSO 1:
        %
        % Calcular o vetor p multiplicando os elementos basicos 
        % de c com a com a matriz inversa das colunas basicas de A
        % p' = cB'*B^(-1) 
        % E tambem calcular os custos reduzidos cj = cj - p'Aj
        % para todos os indices nao-basicos j 
        
	    p = (c(IB(:,1))'*B)';
         
	    j = n+1; % Candidato a entrar na base
	    entrou = 0; % indice no vetor IN que vai entrar na base
        
        printf("Custos reduzidos\n");
        cj = zeros (0, 1);
        
	    for i = 1:rows(IN)
	        cj = [cj; c(IN(i)) - p'*A(:,IN(i)) ];
	        printf("%d %.5f\n", IN(i) , cj(i));
	        
	        % Caso tem algum cj < 0, por Brand vamos
	        % escolher aquele de menor indice
	        if (cj(i) < 0 && IN(i) < j)
	            j = IN(i);
	            entrou = i;
	        endif
	    endfor
	    printf("\n");
        
        % Se j continua valendo n+1 os custos
        % reduzidos sao maior ou igual a zero 
        % a solucao eh % otima, o algortimo para!
        if(j == n+1)
            printf("Solução otima encontrada com custo %.5f\n", c'*x);
            for i = 1:n
                printf("%d %.5f\n", i , x(i));
            endfor
            printf("\n");
            printf("FIM da iteracao ---------------------------\n\n");
            
            ind = 0;
            return;
        endif
        printf("Entra na base: %d\n\n", j);
        
        % Se chegou ate aqui, eh por que tem algum elemento de cj que
        % eh negativo, escolha um deles (usando Bland estamos escolhendo
        % aquele de menor indice, para evitar ciclagem)
        
	    % ------------------------------------------------------------- PASSO 2:
        %
        % Calcula o vetor u = B^(-1)*Aj
	    u = B*A(:,j);
	    
	    printf("Direcao\n");
	    % Ve todos os elementos de u, se nao achar nenhum positivo
	    % o problema eh ilimitado e o problema min c't eh ilimitado
	    ilimitado = 1;
	    for i = 1:rows(u)
	        printf("%d %.5f\n", IB(i) , u(i));
	        if(u(i) > 0)
	            ilimitado = 0;
	        endif
	    endfor
	    printf("\n");
	    
	    % Achamos que o problema eh ilimitado
	    % Mostra a direcao para a qual o problema tem custo -infito
	    if(ilimitado == 1)
	        printf("O problema eh ilimitado na direção da função objetivo\n\n");
	        printf("FIM da iteracao -------------------------------------\n\n");
	        ind = -1;
	        return;
	    endif
	    
	    % ------------------------------------------------------------- PASSO 3:
        %
        % Se existe algum componente u(i) > 0.
        % Define-se:
        % theta* = min[para todo i onde u(i)>0](x(IB(i))/u(i)).
        % k -> um indice tal que theta* = x(IB(k))/u(k).
	    k = 0;
	    theta = 0;
	
	    for i = 1:rows(u)
	        if (u(i) > 0)
	            t = x(IB(i))/u(i);
                if(theta == 0)
                    theta = t;
                    k = i;
                elseif(t < theta)
                    theta = t;
                    k = i;
                endif
	        endif
	    endfor
	    printf("Theta*\n%.5f\n\n", theta);
	    
	    % ------------------------------------------------------------- PASSO 4:
        %
        % Forme uma nova base substituindo A(:,IB(k)) por A(:,j)
        % Os valores da nova solucao basica y sao dados por:
        % y(j) = theta* 
        % y(IB(i)) = x(IB(i)) - theta* u(i), para todo i != k (basicas)
        % y(IB(k)) = 0, (nao-basico)
        
        % Calcula a matriz U, para mudar a nova solucao basica
	    U = zeros (n, 1);
        for i = 1:rows(u)
            U(IB(i)) = -u(i);
        endfor 
        U(j) = 1; % x(j) = 0 + theta*1 = theta*
        
        x = x + theta*U; % Formando um nova solucao basica
        v = x;
               
	    printf("Sai da base: %d\n\n", IB(k));
	    % Coloca IB(k) "aquele que saiu da base", no vetor onde
	    % guarda as variaveis que nao estao na base
	    IN(entrou) = IB(k); 
	    IB(k) = j; % coloca j na base
	    
	    % ------------------------------------------------------------- PASSO 5:
        %
        % Aplicar a sequencia de m operacoes elementares de linhas de 
        % modo a tornar a ultima coluna da matriz [B | u] (que nesse
        % problema nao foi montado, pois nao eh necessaria) igual ao
        % vetor canonico, onde e(k) = 1 e e(i) = 0, para i != k
        
        % Passo 5.1
        % Dividir a k-esima linha de B por u(k), vale recordar que
        % u(k) eh estritamente positiva pois u(k) = -d(IB(k)) > 0.
        % Neste caso estou verificando se u(k) != 1, pois caso seja
        % igual a 1, nao precisamos fazer nenhuma modificacao na
        % linha l da matriz. Economizando m operacoes.
        if(u(k) != 1)
            for i = 1:m
                B(k,i) = B(k,i)/u(k);
            endfor
        endif
        
        % Passo 5.2
        % Para cada linha l de B, para todo l != k, subtrair  
        % k-esima linha de B multiplicada por u(l).
        % Neste caso estou verificando se u(l) != 0, pois caso seja
        % igual a 0, nao precisamos fazer nenhuma modificacao na
        % linha l da matriz. Economizando m operacoes.
        for l = 1:m
            if(l != k && u(l) != 0)
                for i = 1:m
                    B(l,i) = B(l,i) - B(k, i)*u(l);
                endfor
            endif
        endfor
        
        % Fim dessa iteracao
    endwhile
	return;
endfunction
