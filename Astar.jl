#ENV["MPLBACKEND"]="Gtk3Agg"
using DataStructures
#using PyPlot
using ImageView
using Colors
using Gtk







function Extract_M(file)   # Extrait les donnés d'un fichier map pour le convertir en Matrice 
      
    f=open(file) 
    readline(f)
    trans_long::String= split(readline(f)," ")[2]
    trans_larg::String= split(readline(f)," ")[2]
    long=parse(Int, trans_long)
    larg=parse(Int, trans_larg)
    readline(f)
    
    
    
    Mat_ref::Matrix{Char}=Matrix{Char}(undef, long, larg)
    a::Int64=0
    for l in eachline(f)
        a=a+1
       for i = 1:larg
          Mat_ref[a,i]=l[i]
           
       end
    end
    
    return Mat_ref
end 









function poid(v::Vector{Char})         # Retourne le coût de passage entre deux éléments

    function convert_int(x::Char) 
               if (x=='.'||x=='G')
                     return 1
                       elseif (x=='@'||x=='O')
                                 return 2
                           elseif (x=='S')
                              return 3
                           elseif (x=='W')
                              return 4
                                 
                          
               end
    end
    
    c= convert_int(v[1]) 
    d= convert_int(v[2]) 
    Matrice_de_poid::Matrix{Int}=[1 0 1 1;
                                  0 0 0 0;
                                  5 0 5 5;
                                  8 0 8 8]
   
  
    return Matrice_de_poid[d,c]
    
  end 


  

  
  
  function Astar(A::Matrix{Char}, start::Vector{Int64}, final::Vector{Int64} )
  
      
      Visited::Matrix{Int64}= [0 for i in 1:size(A,1), j in 1:size(A,2)] #Matrice : sommet à été visité(true) ou pas visité(false)
      
      Open::Matrix{Int64}= [0 for i in 1:size(A,1), j in 1:size(A,2)] #Matrice : sommet à été visité(true) ou pas visité(false)
  
      Dist_to_start::Matrix{Int64}= [typemax(Int64) for i in 1:size(A,1), j in 1:size(A,2)] #Matrice : Distance de chaque sommet par rapport au start 
  
      Pred::Matrix{Vector{Int64}}= Matrix(undef , size(A,1),size(A,2)) #Matrice : liste des prédécesseurs de chaque sommet
       
      
     
  
  
      function voisins_visitables(som::Vector{Int64}) # Prend les coordonnées d'un sommet et renvoi la liste de ses voisins visitables dans A
        a=[]             
        
         function bornes(a,b) # Fonction auxiliaire qui renvoie vraie si les éléments sont dans la matrice faux sinon
           if a==0 || b==0 || a>=size(A,2)+1 || b>=size(A,1)+1
              return false
             else
               return true
           end
         end
         
         (i,j)=(-1,0)
         k=[som[1]+i,som[2]+j]
         if bornes(k[1],k[2]) && (Visited[k[1],k[2]]==0) && (poid([A[som[1],som[2]],A[k[1],k[2]]])!=0) 
         push!(a,[k[1],k[2]])
         end
  
         (i,j)=(0,-1)
         k=[som[1]+i,som[2]+j]
         if bornes(k[1],k[2]) &&(Visited[k[1],k[2]]==0) && (poid([A[som[1],som[2]],A[k[1],k[2]]])!=0) 
         push!(a,[k[1],k[2]])
         end
  
         (i,j)=(1,0)
         k=[som[1]+i,som[2]+j]
         if bornes(k[1],k[2]) &&(Visited[k[1],k[2]]==0) && (poid([A[som[1],som[2]],A[k[1],k[2]]])!=0) 
         push!(a,[k[1],k[2]])
         end
  
         (i,j)=(0,1)
         k=[som[1]+i,som[2]+j]
         if bornes(k[1],k[2]) &&(Visited[k[1],k[2]]==0) && (poid([A[som[1],som[2]],A[k[1],k[2]]])!=0) 
         push!(a,[k[1],k[2]])
         end
         
        return a
     end 






     function Heuristique(elem,fin)       # Retourne la distance de Manhattan entre deux éléments d'une Matrice
        return abs(fin[1]-elem[1])+abs(fin[2]-elem[2])
     end
  




      # Initialisation 
  
      Visited[start[1],start[2]]=1          
      Dist_to_start[start[1],start[2]]=0
      Pred[start[1],start[2]]=[-1, -1]
      pq= PriorityQueue(start => 0)
      t::Vector{Int64}=[]
      v::Vector{Vector{Int64}}=[]
      cpt=0
  




      
      # Corps de l'algorithme 
         
      while (t!=final) && (!isempty(pq))
        t=dequeue!(pq)
        Visited[t[1],t[2]]=1
        v=voisins_visitables(t)
        cpt=cpt+1
        
        if v==[]
           #on ne fait rien
         else 
          for i in 1:size(v,1)
            
            Open[v[i][1],v[i][2]]=1
            
            if Dist_to_start[t[1],t[2]] + poid([A[t[1],t[2]],A[v[i][1],v[i][2]]])<=Dist_to_start[v[i][1],v[i][2]]
               Dist_to_start[v[i][1],v[i][2]]=Dist_to_start[t[1],t[2]] + poid([A[t[1],t[2]],A[v[i][1],v[i][2]]])
               Pred[v[i][1],v[i][2]]= t
            end
            push!(pq,v[i]=>(Dist_to_start[v[i][1],v[i][2]]+Heuristique(t,final)))
          end 
        end
        
  
      end 
      
      
      



      # Affichage du résultat dans la Matrice 
      
        
        if t==final
          for k in 1:size(A,1)
            for j in 1:size(A,2)
              if Open[k,j]==1
                A[k,j]='Y'
              end

            end
            
          end

          while t!=start
            A[t[1],t[2]]='X'
            t=Pred[t[1],t[2]]
          
          end

         else
          println("Pas de chemin possible")
      end 

    println("longueur du plus court chemin: ",Dist_to_start[final[1],final[2]])
    println("Nombre d'états visités: ",cpt)
  
    return A   
  
end







function char_rgb(c::Char)
    if c == '.'
        return [255, 255, 255]
    elseif c == '@' || c=='O'
        return [0, 0, 0]
    elseif c == 'W'
        return [0, 0, 225]
    elseif c == 'S'
        return [250, 234, 115]
    elseif c == 'X'
        return [255, 0, 0]
    elseif c == 'Y'
        return [121, 128, 129]

    else
        error("Caractère inconnu: $c")
    end
end


function main(nom_fichier,start::Vector{Int64},final::Vector{Int64})



   
    A=Extract_M(nom_fichier)
    
    @time  B=Astar(A,start,final)

  
    
    
    
    color_matrix = fill(RGB(0.,0.,0.),size(B,1),size(B,2))

  
   

    for i in 1:size(B,1)
      for j in 1:size(B,2)
        c = char_rgb(A[i,j])
        color_matrix[i,j]=RGB(c[1]/255.0,c[2]/255.0,c[3]/255.0)
        
      end 
    end 

    
   
    img = imshow(color_matrix)
    resize!(img["gui"]["window"], 800, 800)

    

    
    

end
  