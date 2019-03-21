import numpy as np

class Layer(object):
    def __init__(self, thickness, p_speed, s_speed=None, ps_ratio=None):
        if s_speed is None:
            if ps_ratio is None:
                raise NotImplementedError("Needs either s_speed or ps_ratio")
            s_speed = p_speed / ps_ratio
        self.thickness = thickness
        self.p_speed = p_speed
        self.s_speed = s_speed
        
    def __str__(self):
        return "Layer(thickness={thickness}, p_speed={p_speed}, s_speed={s_speed})".format(
            thickness=self.thickness, p_speed=self.p_speed, s_speed=self.s_speed)
    
    def transmit_p(self, angle_in, x_in, time_in=0):
        hyp = self.thickness / math.cos(math.radians(angle_in))
        time_out = time_in + (hyp / self.p_speed)
        x_out = x_in + (hyp * math.sin(math.radians(angle_in)))
        return x_out, time_out
    
    def transmit_s(self, angle_in, x_in, time_in=0):
        hyp = self.thickness / math.cos(math.radians(angle_in))
        time_out = time_in + (hyp / self.s_speed)
        x_out = x_in + (hyp * math.sin(math.radians(angle_in)))
        return x_out, time_out
    
    def copy(self):
        from copy import deepcopy
        return deepcopy(self)


def snells_law(v1, v2, angle_in):
    return math.degrees(math.asin(v2 * math.sin(math.radians(angle_in))/ v1))


class Model(object):
    def __init__(self, layers):
        """Ordered list of layers, with surface first."""
        if isinstance(layers, list):
            self.layers = layers
        else:
            self.layers = [layers]
        
    def __str__(self):
        return "Model of {0} layers".format(len(self.layers))
    
    def trace_rays(self, source_x, source_depth, angular_interval=0.1):
        layers_above_source = []
        depth = 0.
        for layer in self.layers:
            layers_above_source.append(layer)
            depth += layer.thickness
            if depth >= source_depth:
                bottom_layer = layer.copy()
                bottom_layer.thickness = source_depth - (depth - layer.thickness)
                break
        layers_above_source.reverse()
        
        angular_range = np.arange(
            0, 45. + angular_interval, angular_interval)
        p_times = np.zeros_like(angular_range)
        s_times = np.zeros_like(angular_range)
        p_positions = np.zeros_like(angular_range)
        s_positions = np.zeros_like(angular_range)
        
        for i, angle in enumerate(angular_range):
            x_p, x_s = (source_x, source_x)
            time_p, time_s = (0., 0.)
            angle_p, angle_s = (angle, angle)
            x_p, time_p = bottom_layer.transmit_p(
                angle_in=angle, x_in=source_x, time_in=0)
            x_s, time_s = bottom_layer.transmit_s(
                angle_in=angle, x_in=source_x, time_in=0)
            angle_p = snells_law(
                bottom_layer.p_speed, layers_above_source[1].p_speed, angle)
            angle_s = snells_law(
                bottom_layer.s_speed, layers_above_source[1].s_speed, angle)
            for j, layer in enumerate(layers_above_source[1:]):
                x_p, time_p = layer.transmit_p(
                    angle_in=angle_p, x_in=x_p, time_in=time_p)
                x_s, time_s = layer.transmit_s(
                    angle_in=angle_s, x_in=x_s, time_in=time_s)
                if j < len(layers_above_source) - 2:
                    angle_p = snells_law(
                        layer.p_speed, layers_above_source[j + 2].p_speed, angle_p)
                    angle_s = snells_law(
                        layer.s_speed, layers_above_source[j + 2].s_speed, angle_s)
            p_times[i] = time_p
            s_times[i] = time_s
            p_positions[i] = x_p
            s_positions[i] = x_s
        return angular_range, p_times, p_positions, s_times, s_positions
        
layers = [Layer(thickness=12., p_speed=6.2, ps_ratio=1.68),
          Layer(thickness=11., p_speed=6.6, ps_ratio=1.68),
          Layer(thickness=9.0, p_speed=7.1, ps_ratio=1.68),
          Layer(thickness=19.0, p_speed=8.05, ps_ratio=1.74),
          Layer(thickness=30.0, p_speed=8.25, ps_ratio=1.74),
          Layer(thickness=99999., p_speed=8.5, ps_ratio=1.74)]  # Imaginary deep layer

model = Model(layers)
angular_range, p_times, p_positions, s_times, s_positions = model.trace_rays(
    source_x=0, source_depth=150, angular_interval=1)
fig, axes = plt.subplots(nrows=1, ncols=2)
axes[0].plot(p_times, p_positions, label="P-wave", linestyle="--", marker="x")
axes[1].plot(s_times, s_positions, label="S-wave", linestyle="--", marker="x")
axes[0].set_title("P-wave")
axes[0].set_ylabel("X (km)")
axes[0].set_xlabel("Travel time (s)")
axes[1].set_title("S-wave")
axes[1].set_ylabel("X (km)")
axes[1].set_xlabel("Travel time (s)")