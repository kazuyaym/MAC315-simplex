% EP3 - MAC0315 Programação Linear
% Marcos Kazuya Yamazaki 7577622
% Victor de Queiros Mattoso Archela dos Santos 7557202

%-------------------------------------------------------------------------------------------------
%
%-------------------------------------------------------------------------------------------------

function [ind x d] = simplex(A, b, c, m, n)
    d = [];
    
    printf("***** Simplex: Fase 1 *****\n\n");   
    
    [ind x A IB m] = inicializacaoSimplex(A, b, m, n);
    
    if(ind == 1)
        printf("O problema eh inviavel\n");
        return;
    endif
    
    printf("***** Simplex: Fase 2 *****\n\n"); 
    
    [ind x d] = solucaoSimplex(x, c, A, IB, m, n);
    
    if(ind == 0)
        printf("Solucao otima encontrada com custo: %.3f\n", c'*x);
        x = x'
    else
        printf("O problema eh ilimitado e o custo vai para menos infinito na direcao:\n");
        d = d'
        printf("A partir da solucao viavel:\n");
        x = x'
    endif
	
	return;
endfunction

%-------------------------------------------------------------------------------------------------
%
%-------------------------------------------------------------------------------------------------

% Fase 1 do Simplex: Inicializacao para achar uma solucao viavel
function  [ind xb A IB m] = inicializacaoSimplex(A, b, m, n)
    IB = zeros (0, 1); % indices das variaveis basicas    
    cj = zeros (1, 0);

    % Se necessario, multiplicar as linhas originais por -1, caso b(i) < 0,
    % para obter b >= 0
    for i = 1:m
        if ( b(i) < 0 )
            b(i) = b(i)*(-1);
            for j = 1:n
                A(i,j) = A(i,j)*(-1);
            endfor		    
	    endif
    endfor
    
    % Criar variaveis artificiais y(1:m) 
    I = eye(m);
    A  = [A, I];
    
    % Aplicar o Simplex, com o objetivo de minimizar y(1) + ... + y(m)
    c = [zeros(n,1); ones(m,1)];
    
    for i = 1:m 
        IB = [IB ; n+i];
    endfor
    
    for i = 1:n
        % custos nao basicos, como a matriz basica eh a identidade,
        % os custos valem menos a soma da coluna da matriz A
        cj = [cj, -sum(A(:,i))];
    endfor
    
    cj = [ cj, zeros(1,m) ]; % custos basicos valem zero
    
    [ind x xb A IB d interacao] = metodoSimplex(c, b, cj, A, IB, m, n, 1);
    
    % Se o valor otimo desse problema auxiliar eh positivo, o algoritmo para
    % e o problema original eh inviavel
    if(ind == 0 && c'*x > 0)
        ind = 1;
        return;
    endif
    
    % se o valor otimo for zero, e as variaveis basicas sao todas originais,
    % retire as variaveis auxiliares, e va para a fase 2 
    tirarLinha = zeros(0);
    tirarColuna = zeros(0);
    mudou = false;
    
    for i = 1:m 
        if(IB(i) > n)
            entrou = 0;
            mudou = true;
            
            for j = n:-1:1
                if(A(i,j) != 0)
                    entrou = j;
                endif
            endfor
            
            if(entrou != 0) % coloca 'entrou' na base
                draw(-c'*x, xb, cj, A, IB, m, m+n, n, i, entrou, interacao++, false);
                [A xb] = pivotar(A, xb, i, entrou, m, m+n);
                IB(i) = entrou;
            else % linha LI, retira linha e coluna da base artificial
                tirarLinha = [tirarLinha; i];
                tirarColuna = [tirarColuna; IB(i)];
            endif
        endif
    endfor

    % vamos tirar as colunas e as linhas de tras para frente, 
    % pois caso a matriz tenha 5 colunas, se por exeplo temos que tirar
    % as colunas 1 e 3, se tirarmos a coluna 1 primeiro, a coluna 3 vira
    % a segunda coluna, e assim vamos tirar a coluna 4 sem ter essa intencao.
    tirarColuna = sort(tirarColuna);
    tirarLinha = sort(tirarLinha);
    
    for i = rows(tirarColuna):-1:1;
        draw(-c'*x, xb, cj, A, IB, m, m+n, n, tirarLinha(i), 0, interacao++, false);
        A(:,tirarColuna(i)) = [];
        A(tirarLinha(i),:) = [];
        IB(tirarLinha(i)) = [];
        xb(tirarLinha(i)) = [];
        m = m - 1;
    endfor
    
    % se houve alguma mudanca desde o termino do metodo simplex (retirada de var aux)
    % imprime a tabela nova, caso contrario, nao houve mudancas, e nao eh necessaria 
    % a impressao da tabela
    if(mudou) draw(-c'*x, xb, cj, A, IB, m, m+n, n, 0, 0, interacao, false);
    endif
    
    % Se o valor otimo desse problema auxiliar eh positivo, o algoritmo para
    % e o problema original eh inviavel
    %x = zeros(n+m+rows(tirarColuna), 1);
    %for i = 1:m
    %    x(IB(i)) = xb(i);
    %endfor
    
    
    
	return;
endfunction

%-------------------------------------------------------------------------------------------------
%
%-------------------------------------------------------------------------------------------------

% Fase 2 do Simplex: A partir do ponto obtido pela fase 1,
% acha uma solucao otima, ou nao
function [ind x d] = solucaoSimplex(xb, c, A, IB, m, n);
    % Deletar as variaveis auxiliares e considerar a inversa da base
    for i = (m+n):-1:(n+1)
        A(:,i) = [];
    endfor
    
    % Calcular os custos reduzidos
    cb = zeros(0);
    for i = 1:m
        cb = [cb; c(IB(i))];
    endfor
   
    cj = zeros(0);
    for i = 1:n
        cj = [cj; c(i) - cb'*A(:,i)];
    endfor
    
    % Aplicar o simplex no problema original, a partir da solucao viavel gerada na fase 1
    [ind x xb A IB d interacao] = metodoSimplex(c, xb, cj, A, IB, m, n, 2);
  
	return;
endfunction

%-------------------------------------------------------------------------------------------------
%
%-------------------------------------------------------------------------------------------------

% Execucao do algoritmo do simplex
function [ind x xb A IB d interacao] = metodoSimplex(c, xb, cj, A, IB, m, n, fase)
    d = zeros(0);
    np = n;
    if(fase == 1) n = n+m;
    endif

    interacao = 1;
    while(true)
	    entrou = n+1; % Candidato a entrar na base
        
	    for i = 1:n
	        % Caso tem algum cj < 0, por Brand vamos
	        % escolher aquele de menor indice
	        if (cj(i) < 0 && i < entrou)
	            entrou = i;
	        endif
	    endfor
        
        % Se 'entrou' continua valendo n+1 os custos
        % reduzidos sao maior ou igual a zero 
        % a solucao eh otima, o algortimo para!
        if(entrou == n+1)
            x = zeros(n, 1);
	        for i = 1:m
	            x(IB(i)) = xb(i);
	        endfor
	        
	        draw(-c'*x, xb, cj, A, IB, m, n, np, 0, 0, interacao, true);
            
            ind = 0;
            return;
        endif
        
        %------------------------------------------------------------------
        
        % Ve todos os elementos de u = A(:, entrou), se nao achar nenhum positivo
	    % o problema eh ilimitado e o problema min c'x eh ilimitado
	    
	    % Se existe algum componente A(i, entrou) > 0.
        % Define-se:
        % theta* = min[para todo i onde A(i, entrou) > 0 ](x(IB(i))/A(i, entrou).
        % k -> um indice tal que theta* = x(IB(k))/A(k, entrou).
        
	    k = 0;
	    theta = 0;
	    ilimitado = 1;
	    d = -A(:, entrou);
	    
	    for i = 1:m
	        if(A(i, entrou) > 0)
	            ilimitado = 0;
	            t = xb(i)/A(i, entrou);
                if(theta == 0)
                    theta = t;
                    k = i;
                elseif(t < theta)
                    theta = t;
                    k = i;
                endif
	        endif
	    endfor
	    
	    % Achamos que o problema eh ilimitado
	    % Mostra a direcao para a qual o problema tem custo -infito
	    if(ilimitado == 1)
	        x = zeros(n, 1);
	        for i = 1:m
	            x(IB(i)) = xb(i);
	        endfor
    	    draw(-c'*x, xb, cj, A, IB, m, n, np, 0, 0, interacao, true); 
	        ind = -1;
	        return;
	    endif
        
        %------------------------------------------------------------------
        
        % DESENHAR TABELA!
        x = zeros(n, 1);
        for i = 1:m
            x(IB(i)) = xb(i);
        endfor
        draw(-c'*x, xb, cj, A, IB, m, n, np, k, entrou, interacao, true);
    
        IB(k) = entrou; % coloca 'entrou' na base, e tira o IB(k)
        
        %------------------------------------------------------------------
        
        % Aplicar a sequencia de operacoes elementares de linhas
        [A xb] = pivotar(A, xb, k, entrou, m, n);
        
        % Pivotando a linha dos custos reduzidos
        u = cj(entrou);
        for i = 1:n
            cj(i) = cj(i) - A(k, i)*u;
        endfor
        interacao++;        
    endwhile
	return;
endfunction

%-------------------------------------------------------------------------------------------------
%
%-------------------------------------------------------------------------------------------------

% Fazer pivotacao !
function [A xb] = pivotar(A, xb, k, entrou, m, n)
    % Passo 1
    u = A(k, entrou);
    if(A(k, entrou) != 1)
        for i = 1:n
            A(k, i) = A(k, i)/u;
        endfor
        xb(k) = xb(k)/u; % mexer com o x
    endif

    % Passo 2
    for l = 1:m
        if(l != k && A(l, entrou) != 0)
            u = A(l, entrou);
            for i = 1:n
                A(l,i) = A(l,i) - A(k, i)*u;
            endfor
            xb(l) = xb(l) - xb(k)*u; % mexer com o x 
        endif
    endfor  
    return;
endfunction

%-------------------------------------------------------------------------------------------------
%
%-------------------------------------------------------------------------------------------------

% Funcao para completar com espacos, para nao deixar a tabela deformada 
function drawSpaces (n)
    for i = 1:n
        printf(" ");
    endfor
endfunction

% Funcao para completar com dash
function drawDash (n)
    for i = 1:n
        printf("-");
    endfor
endfunction

function draw(cx, xb, cj, A, IB, m, n, np, lin, col, int, print_custos)
    if(!print_custos && col != 0) printf("Iteracao %d: Forcar a saida de uma variavel auxiliar\n", int);
    elseif(col == 0 && lin != 0) printf("Iteracao %d: Linha redundante, pois a matriz A possui linha(s) LD\n", int);
    else printf("Iteracao %d:\n", int);
    endif

    % transformar -0.000 para 0.000,
    % para que a tabela nao fique deformada
    for i = 1:m
        if(xb(i) == 0) xb(i) = 0;
        endif
        
        for j = 1:n
            if(A(i,j) == 0) A(i,j) = 0;
            endif
        endfor
    endfor
    
    for j = 1:n
        if(cj(j) == 0) cj(j) = 0;
        endif
    endfor

    % Variaveis do problema
    drawSpaces(15);
    printf("|");
    for i = 1:n
        if(i <= np) % numeros de variaveis do problema original
            if(i == col) printf("-> x%d", i);
            else printf("   x%d", i);
            endif
            if(i < 10) drawSpaces(4);
            else drawSpaces(3);
            endif
        else % se n > np, temos que imprimir as variaveis auxiliares
            if(i == col) printf("-> y%d", i-np);
            else printf("   y%d", i-np);
            endif
            if(i-np < 10) drawSpaces(4);
            else drawSpaces(3);
            endif
        endif
        printf("|");
    endfor
    printf("\n");
    
    if(print_custos)
        % Linha zero do tableau, imprimi os custos reduzidos
        if(cx > 99) drawSpaces (7);
        elseif(cx > 9) drawSpaces (8);
        elseif(cx >= 0) drawSpaces (9);
        elseif(cx < -99) drawSpaces (6);
        elseif(cx < -9) drawSpaces (7);
        else drawSpaces (8);
        endif
        
        printf("%.3f |", cx);
      
        for i = 1:n
            if(cj(i) > 99) drawSpaces (1);
            elseif(cj(i) > 9) drawSpaces (2);
            elseif(cj(i) >= 0) drawSpaces (3);
            elseif(cj(i) < -99) drawSpaces (0);
            elseif(cj(i) < -9) drawSpaces (1);
            else drawSpaces (2);
            endif
            printf("%.3f |", cj(i));
        endfor
    else % Caso seja a primeira fase do simplex, na retirada de variaveis auxiliares
         % nao eh preciso mostrar a linha dos custos reduzidos, pois ela nao afeta nas escolhas
        drawSpaces (15);
        printf("|");
        for i = 1:n
            drawSpaces (9);
            printf("|");
        endfor
    endif
    printf("\n");
    
    % dashs
    drawSpaces (3);
    drawDash(12);
    printf("+");
    for i = 1:n
        drawDash(9);
        printf("+");
    endfor
    printf("\n");
    
    % Cada variavel basicas
    for i = 1:m
        if(i == lin) printf("<- ");
        else drawSpaces (3);
        endif
        
        if(IB(i) <= np)
            printf("x%d", IB(i));
            if(IB(i) < 10) drawSpaces (1);
            endif
        else
            printf("y%d", IB(i)-np);
            if(IB(i)-np < 10) drawSpaces (1);
            endif
        endif
        if(xb(i) > 99) drawSpaces (1);
        elseif(xb(i) > 9) drawSpaces (2);
        elseif(xb(i) >= 0) drawSpaces (3);
        elseif(xb(i) < -99) drawSpaces (0);
        elseif(xb(i) < -9) drawSpaces (1);
        else drawSpaces (2);
        endif
        
        printf("%.3f |", xb(i));
        
        for j = 1:n
            if(A(i,j) > 99) drawSpaces (1);
            elseif(A(i,j) > 9) drawSpaces (2);
            elseif(A(i,j) >= 0) drawSpaces (3);
            elseif(A(i,j) < -99) drawSpaces (0);
            elseif(A(i,j) < -9) drawSpaces (1);
            else drawSpaces (2);
            endif
            
            printf("%.3f", A(i,j));
            
            if(i == lin && j == col) 
                printf("*|");
                theta = xb(i)/A(i,j);
            else printf(" |");
            endif
        endfor
        
        printf("\n");
    endfor
    if(lin != 0 && col != 0)
        printf("\nTheta* = %.3f\n\n", theta);
    endif 
    printf("\n");
endfunction

%-------------------------------------------------------------------------------------------------
%
%-------------------------------------------------------------------------------------------------
