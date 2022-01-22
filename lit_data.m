classdef lit_data < handle
  properties
    specs;
    label;
    color;

    % all possible data products that we might get from our literature sources
    Re_s;
    cf;
    alpha;
    G_rat;
    omega;
    mu_torque;
    sigma_torque;
  end
  methods
    function obj = lit_data()

    end
    function init_Ramesh(obj, id)
      fid0 = fopen('./data_directory/torque_data.txt','r');
      C0 = textscan(fid0, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'HeaderLines', 0);
      fclose(fid0);
      ramesh_data = cell2mat(C0);
      ramesh_gap2shear = 2*(17.5)/(17.5 + 16);

      switch id
        case '10'
          % ramesh_G_10_down = ramesh_data(:, 2);
          % ramesh_G_10_up = ramesh_data(:, 5);
          % ramesh_Re_10_down = ramesh_gap2shear*ramesh_data(:, 1);
          % ramesh_Re_10_up = ramesh_gap2shear*ramesh_data(:, 4);
          % ramesh_Grat_10_down = ramesh_data(:, 8);
          ramesh_Grat_10_up = ramesh_data(:, 11);
          % ramesh_Re_4Grat_10_down = ramesh_gap2shear*ramesh_data(:, 7);
          ramesh_Re_4Grat_10_up = ramesh_gap2shear*ramesh_data(:, 10);

          obj.Re_s = ramesh_Re_4Grat_10_up;
          obj.G_rat = ramesh_Grat_10_up;
          obj.label = 'RM10';
          obj.specs = ' -';
          obj.color = [255, 191, 0]/255;
        case '20'
          % ramesh_G_20_down = ramesh_data(:, 14);
          % ramesh_G_20_up = ramesh_data(:, 17);
          % ramesh_Re_20_down = ramesh_gap2shear*ramesh_data(:, 13);
          % ramesh_Re_20_up = ramesh_gap2shear*ramesh_data(:, 16);
          % ramesh_Grat_20_down = ramesh_data(:, 20);
          ramesh_Grat_20_up = ramesh_data(:, 23);
          % ramesh_Re_4Grat_20_down = ramesh_gap2shear*ramesh_data(:, 19);
          ramesh_Re_4Grat_20_up = ramesh_gap2shear*ramesh_data(:, 22);

          obj.Re_s = ramesh_Re_4Grat_20_up;
          obj.G_rat = ramesh_Grat_20_up;
          obj.label = 'RM20';
          obj.specs = ' -';
          obj.color = [191, 0, 255]/255;
      end
    end
    function init_Racina(obj)
      run figure_properties.m
      re_racina_a = 800:10:1e4 ;
      cf_racina_a = Racina_gen_cf(re_racina_a, 'a');
      % re_racina_b = 1e4:100:1e6 ;
      % cf_racina_b = Racina_gen_cf(re_racina_b, 'b');
      % log_re_racina_a = log(re_racina_a);
      % log_cf_racina_a = log(cf_racina_a);
      % log_cf_racina_b = log(cf_racina_b);
      % primelogRK_cf_re_racina_a = approx_deriv_2ndO_legrangian(log_re_racina_a, log_cf_racina_a);
      % alphaRK_racina_a = primelogRK_cf_re_racina_a + 2;

      obj.Re_s = re_racina_a;
      obj.cf = cf_racina_a;
      obj.label = 'RK';
      obj.specs = ':';
      obj.color = purple1;

    end
    function init_Ravelet(obj, id)
      obj.color = [0 0 0];
      obj.specs = ' .';
      obj.label = 'RV';
      switch id
        case 'alpha' % lower range of Reynolds numbers
          data_Ravelet = dlmread('./data_directory/DataRavelet_alphaFig5in.txt');
          obj.Re_s = data_Ravelet(:, 1);
          obj.alpha = data_Ravelet(:, 2);
        case 'cf' % higher range of Reynolds numbers
          data_Ravelet = dlmread('./data_directory/DataRaveletFig5in.txt');
          obj.Re_s = data_Ravelet(:, 1);
          obj.cf = data_Ravelet(:, 2);
      end
    end
    function init_Lewis(obj, id)
      run figure_properties.m
      %% conversion factor for gap Reynolds number to shear Reynolds number in the case of Lewis and Swinney
      gap2shear_Lewis = 2*(22.085e-2)/( ( 22.085e-2 ) + ( 15.999e-2 ) );

      %% Lewis and Swinney produce two polynomial fits, valid for two separate ranges of gap Reynolds numbers
      switch id
        case 'a' % lower range of Reynolds numbers
          Reb = ((2000:1000:11000));
          log10_Reb = log(Reb)/log(10);

          obj.alpha = 3*(0.2005)*(log10_Reb).^(2) + 2*(-1.970)*(log10_Reb) + 7.775;
          obj.label = 'LSa';
          obj.specs = ' -';
        case 'b' % higher range of Reynolds numbers
          Reb = ((13000:1000:1e6));
          log10_Reb = log(Reb)/log(10);

          obj.alpha = 3*(-0.006360)*(log10_Reb).^(2) + 2*(0.1349)*(log10_Reb) + 0.8850;
          obj.label = 'LSb';
          obj.specs = ' --';
      end
      obj.Re_s = gap2shear_Lewis*Reb;
      obj.color = orange4;
    end
    function init_Chuan(obj, id)
      run figure_properties.m
      fluidization_colours4 = [255, 56, 214; 242, 52, 207; 229, 48, 200; 217, 44, 194; 204, 40, 187; 192, 36, 181; 179, 32, 174; 167, 28, 168; 154, 24, 161; 141, 20, 154; 129, 16, 148; 116, 12, 141; 104, 8, 135; 91, 4, 128; 79, 0, 122]*(255)^(-1);
      obj.specs = 'h';
      switch id
        case '0.0'
          obj.omega = [523.596787544967e-003; 10.4719317526646e-003; 104.718430589292e-003; 5.23555843046626e+000; 10.4710786569978e+000; 52.3555698229099e+000; 104.715448666809e+000];
          obj.mu_torque = [6.98771449515448e-003; 7.69096692200915e-003; 7.60614292222222e-003; 7.50320827250608e-003; 7.24839593679458e-003; 7.61510178173719e-003; 8.09913595505618e-003];
          obj.sigma_torque = [53.4906393458386e-006; 408.314499239035e-006; 70.3662817579475e-006; 34.6602665561042e-006; 27.2589996217045e-006; 110.616728288155e-006; 71.1017514881501e-006];
          obj.label = 'FB1 q=0.0';
          obj.color = fluidization_colours4(1, :);
        case '0.5'
          obj.omega = [10.4709010659071e+000; 20.9399438904452e+000; 31.4116957049155e+000; 41.8818824189874e+000; 52.3530884648949e+000; 62.8615129285921e+000; 73.2961174017939e+000; 83.7691695018464e+000; 94.2413796740110e+000; 104.713836647983e+000; 115.184393124975e+000];
          obj.mu_torque = [6.28525765306123e-003; 6.28879863945578e-003; 6.28362700729927e-003; 6.56845939597315e-003; 6.79788006756757e-003; 6.98355538461538e-003; 7.34339259259259e-003; 7.72254194630873e-003; 8.13137466216216e-003; 8.60804931034482e-003; 9.03217517006803e-003];
          obj.sigma_torque = [44.5806695214019e-006; 38.0139034855783e-006; 35.8602875102105e-006; 115.299408782644e-006; 147.368858937030e-006; 25.3092055545984e-006; 133.769450766827e-006; 118.865955198661e-006; 127.523995747510e-006; 71.0284042033536e-006; 36.2119061342456e-006];
          obj.label = 'FB1 q=0.3';
          obj.color = fluidization_colours4(5, :);
        case '0.75'
          obj.omega = [10.4690873336664e+000; 20.9398989811737e+000; 31.4103473786741e+000; 41.8836885762354e+000; 52.3542572470467e+000; 62.7937343757488e+000; 73.3008562351569e+000; 83.7718621197607e+000; 94.2419811050370e+000; 104.714018606973e+000; 115.182473834160e+000];
          obj.mu_torque = [4.47530925266904e-003; 4.54697869415808e-003; 4.69721993243243e-003; 5.00056498316498e-003; 5.29476700336700e-003; 5.85891784511785e-003; 6.17057070707071e-003; 6.70186599326599e-003; 7.29465844594595e-003; 7.92223661016949e-003; 8.45692662116041e-003];
          obj.sigma_torque = [42.8962425903257e-006; 34.9826956190008e-006; 93.2249503301729e-006; 112.883472203444e-006; 135.927976216652e-006; 35.2420776712163e-006; 123.581684540591e-006; 108.496938967469e-006; 96.7370046260344e-006; 63.4278958729085e-006; 34.1307815291551e-006];
          obj.label = 'FB1 q=0.4';
          obj.color = fluidization_colours4(9, :);
        case '1.0'
          obj.omega = [10.4699789320195e+000; 20.9395442164895e+000; 31.4124990189983e+000; 41.8819195780697e+000; 52.3541108976392e+000; 62.8128732194048e+000; 73.2969934970426e+000; 83.7704938135510e+000; 94.2413730152678e+000; 104.715806901865e+000; 115.181944751061e+000];
          obj.mu_torque = [2.60611904761905e-003; 2.74518225255973e-003; 2.97735904436860e-003; 3.34200135135135e-003; 3.74963581081081e-003; 4.16402764976959e-003; 4.90227027027027e-003; 5.62718817567568e-003; 6.38328181818182e-003; 7.10244560810811e-003; 7.73245310344828e-003];
          obj.sigma_torque = [32.8655662461165e-006; 35.3587750299827e-006; 11.4620888705079e-006; 105.482330104475e-006; 118.406867935334e-006; 10.3007492763657e-006; 113.043361687547e-006; 99.9893033027316e-006; 104.721787496264e-006; 58.0413346074914e-006; 27.5369950721156e-006];
          obj.label = 'FB1 q=0.5';
          obj.color = fluidization_colours4(13, :);
      end
    end
  end
end

function prime_vec = approx_deriv_2ndO_legrangian(t, x)
    prime_vec = zeros(1, length(t)-4 );
    for i = 1:length(prime_vec)
        c = polyfit(t(i:i+3),x(i:i+3),1);
        prime_vec(i) = c(1);
    end
end

function cf = Racina_gen_cf( Res, tag )
  % so, when we use our eta for both, we have our best fit. When we use
  % their eta for shear2gap, it may be more accurate. not sure
  eta = 0.4832; % our eta
  shear2gap = 0.5*(1 + eta); % so, the idea is as below:
  % Racina gave us a way of plotting approximate G for Reb. We take OUR
  % Reb, which is calculated by OUR eta and OUR Res, hand it off to their
  % model relating expected G for a given Reb, then turn that into a
  % coefficient of friction. So I think this is right.
  Reb_vec = Res*shear2gap;
  cf = zeros(1, length(Res) );
  C = (2*pi)^(-1)*( 1 - eta^2 )^(2)*(2*eta)^(-2);
  for i = 1:length(Res)
      Reb = Reb_vec(i);
      prod1 = (Reb)^(-2);
      if tag == 'a'
          prod2 = 2.13 * (eta)^(3/2) * (1-eta)^(-7/4)*(Reb)^(1.445);
      elseif tag == 'b'
          prod2 = 0.113 * (eta)^(3/2) * (1-eta)^(-7/4)*(Reb)^(1.764);
      end
      cf(i) = C*prod1*prod2;
  end
end
