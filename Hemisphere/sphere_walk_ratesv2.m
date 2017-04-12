function  sphere_walk_ratesv2(seconds,tsize,cells,mode)

%N.B. This is different to gravity_vector, as no longer saving every 100th
%timestep....

    input_time = uint32(((0.1/tsize)));
    steps = uint32(((seconds/tsize)));
    path_time = steps-input_time;
 
        
    x = (tsize:(tsize):seconds); %Now creates x as going from tsize to seconds in steps of tsize.
    xtick = [0,steps/4, steps/2, (steps/4)*3, steps];
    xticklabel = [0,seconds/4, seconds/2, (seconds/4)*3, seconds];
    


y1 = (1:cells);
[X1, Y1] = meshgrid(x, y1);

         
%LOADING CELL FIRING RATES%             
HDRates = zeros(cells(1), steps);
HDRates = file_load(cells(1), steps, 'HDRates.bdat');

yaw_input = zeros(cells(1), steps);
yaw_input = file_load(cells(1), steps, 'combined_yaw_input.bdat');

HD_pvector = zeros(1,path_time);
HD_pvector = file_load(1,path_time,'Pvector.bdat');

Input_pvector = zeros(1,path_time);
Input_pvector = file_load(1,path_time,'InputPvector.bdat');

                
%PLOTTING CELL RATES AND INPUTS%
figure();
grid on;

subplot(2,1,1);

    h1 = surf(X1, Y1, HDRates);

    set(h1, 'LineStyle', 'none');
    xlabel('Time (s)', 'Fontsize', 10);
    %set(gca, 'XTick', xtick, 'Fontsize', 24);
    %set(gca, 'XTickLabel',xticklabel,'Fontsize',24);
    ylabel('Cell', 'Fontsize', 10);
    ylim([1 cells(1)]);
    %xlim([0 seconds]);
    colormap(1-gray);
    view(0, 90);
    title('HD Cell Firing');
    
    subplot(2,1,2);

    h1 = surf(X1, Y1, yaw_input);

    set(h1, 'LineStyle', 'none');
    xlabel('Time (s)', 'Fontsize', 10);
    %set(gca, 'XTick', xtick, 'Fontsize', 24);
    %set(gca, 'XTickLabel',xticklabel,'Fontsize',24);
    ylabel('Cell', 'Fontsize', 10);
    ylim([1 cells(1)]);
    %xlim([0 seconds]);
    colormap(1-gray);
    view(0, 90);
    title('Yaw Input');
    
    if (strcmpi(mode, 'save'))
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    saveas(gcf,'_rates', 'fig');
    close(gcf);
    end
      
    
   
    figure()
    plot(HD_pvector(1:1000),'b');
    hold on
    plot(Input_pvector(1:1000),'r');
    xlabel('Path Timestep', 'Fontsize', 10);
    ylabel('Position', 'Fontsize', 10);
    %xlim([0 seconds]);
    ylim([0 360]);
    title('HD and Input locations');
    legend('HD Location', 'Input Location');
  
    
    
             
                
end


function rates = file_load(cells, steps, fname)
        
       rates = zeros(cells, steps);
       
       fid = fopen(fname, 'rb');
       
       rates = fread(fid, [steps, cells], 'float32')';
       
       fclose(fid);
       
  
end


