%% CURRENT LATEST VERSION (01/10/23)

clc;
clear;

while true
n_blocks = input('Enter the number of components: ');
  if isempty(n_blocks)
    disp('Error: You must enter the number of components')
  else
   break;
  end
end
%%inputting lengths, weights, stiffnesses, CGs

for i = 1:n_blocks

    %defining block properties and making sure all the inputs are taken in
    while true
    block_length = input(['Input the length of the component ' num2str(i) ': ']);
      if isempty(block_length)
        disp('Error: You must enter the component length.');  % Display an error message if the input is empty
      else
       block_length(i) = block_length;
       break;
      end
     end

    while true
    block_mass = input(['Enter the mass of component ' num2str(i) ': ']);
      if isempty(block_mass)
        disp('Error: You must the component mass.');  % Display an error message if the input is empty
      else
        block_mass(i) = block_mass;
        break;
      end
    end

    while true
    block_cg = input(['Enter the distance of the center of gravity for component ' num2str(i) ' from the output flange moment in m: ']);
      if isempty(block_cg)
          disp('Error: You must enter the component center of gravity.');  % Display an error message if the input is empty
      else
        block_cg(i) = block_cg;
        break;
      end
    end

end

cumul = 0;

for i=1:length(block_length)
  cumul = cumul + block_length(i);
  cumul_block(i) = cumul;
endfor

%%inputting spring properties

while true
n_mounts = input('Enter the number of supports: ');
  if isempty(n_mounts)
    disp('Error: You must enter the number of supports')
  else
   break;
  end
end


for i = 1:n_mounts

  %defining mount and bracket properties
  while true
    mount_distance = input(['Input the distance of support ' num2str(i) ' from the output flange moment in m: ']);
      if isempty(mount_distance)
        disp('Error: You must enter the support distance.');  % Display an error message if the input is empty
      else
       mount_distance(i) = mount_distance;
       break;
      end
  end

  while true
    bracket_stiffness = input(['Input the stiffness of bracket ' num2str(i) ' in N/m:']);
      if isempty(bracket_stiffness)
        disp('Error: You must enter the bracket stiffness.');  % Display an error message if the input is empty
      else
       bracket_stiffness(i) = bracket_stiffness;
       break;
      end
  end

  while true
    mount_stiffness = input(['Input the stiffness of the mount used for bracket ' num2str(i) ':']);
      if isempty(mount_stiffness)
        disp('Error: You must enter the mount_stiffness.');  % Display an error message if the input is empty
      else
       mount_stiffness(i) =  mount_stiffness;
       break;
      end
  end

  while true
     mount_bottoming_out = input(['Input the bottoming out displacement for the mount used for bracket ' num2str(i) ':']); %displacement till bottomout
      if isempty(mount_bottoming_out)
        disp('Error: You must enter the bottoming-out displacement for the mount.');  % Display an error message if the input is empty
      else
       mount_bottoming_out(i) = mount_bottoming_out;
       break;
      end
   end

  support_stiffness(i) = 2* (mount_stiffness(i)*bracket_stiffness(i)/(mount_stiffness(i)+bracket_stiffness(i)));
  bottoming_out(i) =  mount_bottoming_out(i)*(1+(mount_stiffness(i)/bracket_stiffness(i)));

 end

 og_bottoming_out = bottoming_out;
%% output flange moment and g-force

while true
output_flange_moment = input('Enter the moment at the output flange: ');
  if isempty(output_flange_moment)
    disp('Error: You must enter the moment at the output flange. If there is none, input 0.')
  else
   break;
  end
end

while true
g_s = input('Enter the acceleration experienced by the system in m/(s^2): ');
  if isempty(g_s)
    disp('Error: You must enter the g-force value.')
  else
   break;
  end
end


%% calculating the effective center of gravity from first spring and total mass
net_cg = sum(block_mass.*block_cg)/sum(block_mass)-mount_distance(1);
total_mass = sum(block_mass);
mass_moment = total_mass*net_cg;
L = mount_distance(n_mounts) - mount_distance(1);

%% calculating stiffness matrix components

function stiffness_matrix = calculateStiffnessMatrix(support_stiffness, mount_distance, L, n_mounts)

k11 = support_stiffness(1);
k12 = support_stiffness(n_mounts);
k21 = 0;
k22 = support_stiffness(n_mounts)*L;

for i = 1:n_mounts-2

  k_add11 = support_stiffness(i+1)*(1-(mount_distance(i+1)-mount_distance(1))/L);
  k11 = k11 + k_add11;

  k_add12 = support_stiffness(i+1)*(mount_distance(i+1)-mount_distance(1))/L;
  k12 = k12 + k_add12;

  k_add21 = support_stiffness(i+1)*(mount_distance(i+1)-mount_distance(1))*(1-(mount_distance(i+1)-mount_distance(1))/L);
  k21 = k21 + k_add21;

  k_add22 = support_stiffness(i+1)*(mount_distance(i+1)-mount_distance(1))*(mount_distance(i+1)-mount_distance(1))/L;
  k22 = k22 + k_add22;

 end
stiffness_matrix = [k11,k12;k21,k22];
end

%assembling force matrices
mass_matrix = [total_mass;mass_moment];
external_force_matrix = [0;output_flange_moment];
variable_weight_matrix = g_s.*mass_matrix;

%%solving the initial problem
%making the pilot stiffness matrix
pilot_stiffness_matrix = calculateStiffnessMatrix(support_stiffness,mount_distance,L,n_mounts);

%external moment displacement
output_moment_displacement = inv(pilot_stiffness_matrix)*(external_force_matrix);

%self weight displacement
self_weight_displacement = inv(pilot_stiffness_matrix)*(variable_weight_matrix);

%calculating displacements in springs at the ends
pilot_displacement_matrix = self_weight_displacement + output_moment_displacement;

%calculating displacements in springs other than the first and last springs
 for i = 1:n_mounts-2

   pilot_disp(i) = pilot_displacement_matrix(1)+(pilot_displacement_matrix(2)-pilot_displacement_matrix(1))*mount_distance(i+1)/L;

 end

 %making full displacement vector
pilot_displacements = horzcat(pilot_displacement_matrix(1),pilot_disp,pilot_displacement_matrix(2));

 %checking if there are bottom outs
difference = bottoming_out - pilot_displacements;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NO BOTTOM OUT CASE SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(difference <= 0) == 0
  total_displacement = pilot_displacements;
  disp('Solution reached - displacements:');
  disp(total_displacement);
  reaction_force = support_stiffness.*total_displacement;
  disp('Solution reached - reaction forces:')
  disp(reaction_force)
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BOTTOM OUT CASE SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(difference <= 0) == 1 %if there are bottom-outs, this is the next step

delta_g = [];
del_g = g_s;                             %this is the value that the load starts at and is reduced every loop until there are no bottom outs
% Define symbolic variables                          %the symbolic variables defines here are displacement values for the first and last supports and the g-force acting on the system
stiffness_matrix = pilot_stiffness_matrix; %initialises the stiffness matrix with the pilot stiffness matrix from earlier. this will be changes as the iterations continue
flag = 0;                                  %sets a flag as false that is used to separate the first bottom out load from the subsequent loads
displacement_beyond_current_BO = 0;          %initialises
true_displacement_at_1st_BO = 0;           %initializing the real result, which will be added to it iteratively
BO_checker = zeros(1,n_mounts);
true_displacement=0;
stiffness_calc = support_stiffness;

itrnum = 0; %starting value for iteration number

qqq = [];

  while any(difference <= 0) == 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MINIMUM G-FORCE MOUNT SOLVING MODULE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    itrnum = itrnum +1;
disp(['ITERATION ', num2str(itrnum)]); %keeps track of iterations

 % Solve the equations for the two cases where displacements at the ends are known
    % Case where x1 is fixed
    A1 = [stiffness_matrix(1,2), -(mass_matrix(1));stiffness_matrix(2,2), -(mass_matrix(2))];
    b1 = [(-stiffness_matrix(1,1)*bottoming_out(1)+external_force_matrix(1));(-stiffness_matrix(2,1)*bottoming_out(1)+external_force_matrix(2))];
    x1 = A1\b1;
    g_sol1 = x1(2);


    % Case where x2 is fixed
    A2 = [stiffness_matrix(1,1), -(mass_matrix(1));stiffness_matrix(2,1), -(mass_matrix(2))];
    b2 = [(-stiffness_matrix(1,2)*bottoming_out(n_mounts)+external_force_matrix(1));(-stiffness_matrix(2,2)*bottoming_out(n_mounts)+external_force_matrix(2))];
    x2 = A2\b2;
    g_soln = x2(2);


  % Solve the equations for x1, x2, and g in each case

      % Define the known values of z
      z = bottoming_out(2:end-1); % Add more values if needed
      for i = 1:length(z)
          %B_x = [(z(i)-x2*sym(mount_distance(i+1)/L))/(1-sym(mount_distance(i+1)/L));((z(i)-x1)*sym(L/mount_distance(i+1))+x1)]
          A3 = [stiffness_matrix(1,1),stiffness_matrix(1,2), (-mass_matrix(1));stiffness_matrix(2,1), stiffness_matrix(2,2), (-mass_matrix(2)); (1-mount_distance(i+1)/L), (mount_distance(i+1)/L),0];
          b3 = [external_force_matrix(1);external_force_matrix(2);z(i)];
          x3 = A3\b3;

          g_sol(i) = x3(3);
      endfor

    g_bottoming_out = horzcat(g_sol1,g_sol,g_soln); %finds what causes g each mount to bottom out for the current case

    for i = 1:length(BO_checker)
      if BO_checker(i) == 1
        g_bottoming_out(i) = inf;
      endif
    endfor

    positive_indices = find(g_bottoming_out > 0);                                                 % This section takes the g_bottoming_out array
    [first_bottoming_out_g, first_bottoming_out_mount] = min(g_bottoming_out(positive_indices));  % and extracts the minimum positive value for
    first_bottoming_out_mount = positive_indices(first_bottoming_out_mount);                   % a mount to bottom out and its index % If it's negative, it means that the force is upwards

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATING DISAPLACEMENTS AT FIRST BOTTOM OUT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Note: 1st bottoming out here refers to the first mount to bottom out in the context of this iteration, not in terms of the whole problem.

    displacements_1st_BO = inv(stiffness_matrix)*(first_bottoming_out_g.*mass_matrix) + output_moment_displacement; %calculates displacements upto the bottom out of the given spring
    %calculating displacements in springs other than the first and last springs

    for i = 1:(n_mounts-2)
      disp_1stBO(i) = displacements_1st_BO(1)+(displacements_1st_BO(2)-displacements_1st_BO(1))*mount_distance(i+1)/L; %calculates the displacements of the mounts in between the first and last ones
    endfor
    %making full displacement vector
    true_displacement_at_1st_BO  = horzcat(displacements_1st_BO(1),disp_1stBO,displacements_1st_BO(2));

    if flag == 0
    disp(['The first mount bottoms out at ', num2str(first_bottoming_out_g), ' N/m2']); %Extracting g force for the very first bottom out, considered salient
    flag = 1;
  endif

%% Note : The above 'if' statement uses a flag to ensure that the value of first_bottoming_out_g recorded is only the very first one and none other

    true_displacement = true_displacement + true_displacement_at_1st_BO; %iteratively adding the true displacement at current BO level

%%%%%%%%%%%%%%%%%%%%%%%%%  RECALCULATING THE DISPLACEMENTS BEYOND THE CURRENT BOTTOM OUT TO CHECK FOR FURTHER BOTTOM OUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%


    BO_checker(first_bottoming_out_mount) = 1; %updates the check array to say that the given mount has bottomed out

    stiffness_calc(first_bottoming_out_mount) = 2*bracket_stiffness(first_bottoming_out_mount); %switches out the stiffness of the support that bottoms out to only that of the bracket at that point    new_stiffness_matrix = calculateStiffnessMatrix(support_stiffness, mount_distance, L, n_mounts);
    stiffness_matrix = calculateStiffnessMatrix(stiffness_calc,mount_distance,L,n_mounts); %calculating the stiffness matrix


    del_g = del_g - first_bottoming_out_g; %finding the remaining load to be applied after the current bottom out. this decreases as the amount of load to be applied reduces. the final del_g is the one where there are no more bottom outs.
    qqq = [qqq,del_g];

    delta_g = [delta_g,first_bottoming_out_g]; %creates an array that keeps track of all the deltas beyond the first bottom out for subsequent bottom outs

    extra_displacement= inv(stiffness_matrix)*(del_g.*mass_matrix); %calculates the displacement beyond the first bottom out. This could cause bottom outs.


    for i = 1:n_mounts-2
      rem_check_disp(i) = extra_displacement(1)+(extra_displacement(2)-extra_displacement(1))*mount_distance(i+1)/L; %finds the values of displacements between the mounts at the ends
    endfor
    displacement_beyond_BO  = horzcat(extra_displacement(1),rem_check_disp,extra_displacement(2));

    total_displacement = true_displacement+displacement_beyond_BO; %applies all of the remaining g force and it used to check for further bottom outs
    disp_checker = true_displacement_at_1st_BO+displacement_beyond_BO;
    for i=1:length(BO_checker)
      if BO_checker(i) == 1
        disp_checker(i) = -inf;
      endif
    endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECKING RESULTS FOR FURTHER BOTTOM OUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    difference = bottoming_out - disp_checker; %checking if there are any more bottom outs

    if all(difference >=0) == 1
      disp('Solution reached - displacements');
      disp(total_displacement);
      break                     %If all differences are greater than or equal to zero, this means there are no more mounts bottoming out and the final solution has been reached
    endif

    for i=1:length(BO_checker)  %for all the mounts that haven't bottomed out, we check how much further they can go before they bottom out
      if BO_checker(i) == 0
        bottoming_out(i) = bottoming_out(i) - true_displacement_at_1st_BO(i);
      endif
    endfor

endwhile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FINDING REACTION FORCES PER SUPPORT IN B/O CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:n_mounts %done in a loop to take each support one by one and calculate the relevant reaction force and put them into an array
  if (total_displacement(i)>og_bottoming_out(i)) == true
    reaction_force(i) = support_stiffness(i)*og_bottoming_out(i)+2*bracket_stiffness(i)*(total_displacement(i)-og_bottoming_out(i)); %calculates reaction forces by splitting the displacement into a pre and post bottom out scenario and applying the relevant stiffnesses
  endif
  if (total_displacement(i)>og_bottoming_out(i)) == false
    reaction_force(i) = support_stiffness(i)*total_displacement(i); %calculates the reaction forces by simply multiplying the displacement by support stiffness
  endif
endfor
disp('Solution reached - reaction forces:')
disp(reaction_force/2)
endif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FINDING BENDING MOMENTS AT THE FLANGES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:n_mounts %done in a loop to take each support one by one and calculate the relevant reaction force and put them into an array
  if (total_displacement(i)>og_bottoming_out(i)) == true
    reaction_force(i) = (support_stiffness(i)*og_bottoming_out(i)+2*bracket_stiffness(i)*(total_displacement(i)-og_bottoming_out(i))); %calculates reaction forces by splitting the displacement into a pre and post bottom out scenario and applying the relevant stiffnesses
  endif
  if (total_displacement(i)>og_bottoming_out(i)) == false
    reaction_force(i) = (support_stiffness(i)*total_displacement(i)); %calculates the reaction forces by simply multiplying the displacement by support stiffness
  endif
endfor
disp('Solution reached - reaction forces:')
disp(reaction_force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FINDING BENDING MOMENTS AT THE FLANGES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (length(block_length)>1) == true
itr=0;
for i=1:(length(block_length)-1) %makes a loop for each flange
  itr=itr+1;
  
  moment_left=0;
  moment_right=0;
  
  for j =1:length(block_length)  %loops through all the blocks and sums up moments on either side of the flange from the blocks
    if block_cg(j) < cumul_block(i) == true
      moment_left = moment_left+block_mass(j)*g_s*abs(block_cg(j)-cumul_block(i)); %adds moments on the left of the flange caused by the block weights
    endif

    if block_cg(j) < cumul_block(i) == false
      moment_right = moment_right-block_mass(j)*g_s*abs(block_cg(j) - cumul_block(i));  %adds moments on the right of the flange caused by the block weights
    endif
  endfor

  for kk=1:length(mount_distance)
        if mount_distance(kk) < cumul_block(i) == true
          moment_left = moment_left-reaction_force(kk)*abs(mount_distance(kk)-cumul_block(i));  %adds moments on the left of the flange caused by the support reactions
        endif
        if mount_distance(kk) > cumul_block(i) == true
          moment_right = moment_right + reaction_force(kk)*abs(mount_distance(kk)-cumul_block(i)); %adds moments on the right of the flange caused by the support reactions
        endif
  endfor
 flange_moments(i) = moment_left - moment_right; %takes the resultant moment from the two sides. it's subtraction here because different sign conventions were used on either side
endfor

disp('The moments at the flanges respectively are:');
disp(flange_moments); %displays the moments at the flanges
endif

if (length(block_length)>1) == false
  disp('No flange moments between components are there is only one component') %in case there is only one component
endif
