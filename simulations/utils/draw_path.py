import numpy as np
import matplotlib.pyplot as plt
import pickle
import scipy.io
import math
import copy
import cv2
import time
from threading import Thread, Lock
from typing import Union, Tuple, List
import sys

class MovementPrimitive:

    def __init__(self, n_dofs: int, n_kernels: int, kernels_std_scaling: float = 1.5):
        self.__init_helper__(n_dofs, n_kernels, kernels_std_scaling)

    def __init_helper__(self, n_dofs: int, n_kernels: int, kernels_std_scaling: float = 1.5):
        self._config = {'n_dofs': n_dofs, 'n_kernels': n_kernels, 'kernels_std_scaling': kernels_std_scaling}
        self.n_dofs = n_dofs
        self.n_kernels = n_kernels
        self.kernels_std_scaling = kernels_std_scaling

        self.weights = np.zeros((n_dofs, n_kernels))

        self.c = np.linspace(0.0, 1.0, self.n_kernels)[:, np.newaxis]

        hi = 1.0 / (kernels_std_scaling * (self.c[1] - self.c[0])) ** 2
        self.h = np.ones((n_kernels, 1)) * hi

        zero_tol = np.nextafter(np.float32(0), np.float32(1))
        self.s_min = self.c[0] - math.sqrt(-math.log(zero_tol) / self.h[0])
        self.s_max = self.c[-1] + math.sqrt(-math.log(zero_tol) / self.h[-1])

    def get_pos(self, s: float) -> np.array:
        return np.matmul(self.weights, self.regress_vec(s))

    def get_vel(self, s: float, s_dot: float) -> np.array:
        return np.matmul(self.weights, self.regress_vec_dot(s, s_dot))

    def get_accel(self, s: float, s_dot: float, s_ddot: float) -> np.array:
        return np.matmul(self.weights, self.regress_vec_ddot(s, s_dot, s_ddot))

    def train(self, s: np.array, pos_data: np.array, train_method: str = 'LS', end_points_constraints=False):

        s = np.squeeze(s)

        if pos_data.ndim == 1:
            pos_data = np.expand_dims(pos_data, 0)

        if pos_data.shape[0] != self.n_dofs:
            raise AttributeError('[MovementPrimitive::train]: The training data have wrong number of DoFs...')

        if np.any(s > 1) or np.any(s < 0):
            print('\33[1m\33[33m[MovementPrimitive::train]: The training timestamps are not normalized...\33[0m')

        H = np.hstack([self.regress_vec(s[j]) for j in range(len(s))])

        if train_method.upper() == 'LWR':
            self.weights = np.matmul(pos_data, H.transpose()) / np.sum(H, axis=1).transpose()
        elif train_method.upper() == 'LS':
            self.weights = np.linalg.lstsq(H.transpose(), pos_data.transpose(), rcond=None)[0].transpose()
        else:
            print('\33[1m\33[31m[MovementPrimitive::train]: Unsupported training method...\33[0m')

        if end_points_constraints:
            Sw = np.linalg.inv(np.matmul(H, H.transpose()))

            # enforce start and end point constraints
            A = np.hstack([self.regress_vec(0), self.regress_vec(1),
                           self.regress_vec_dot(0, 1), self.regress_vec_dot(1, 1),
                           self.regress_vec_ddot(0, 1, 0), self.regress_vec_ddot(1, 1, 0)])
            b = np.stack([pos_data[:, 0], pos_data[:, -1],
                          np.zeros((self.n_dofs,)), np.zeros((self.n_dofs,)),
                          np.zeros((self.n_dofs,)), np.zeros((self.n_dofs,))],
                         axis=1)

            R = np.diag([1e-5, 1e-5, 1e-3, 1e-3, 1e-2, 1e-2])
            Sw_At = np.matmul(Sw, A)  # A is already transposed!
            K = np.matmul(Sw_At, np.linalg.inv(R + np.matmul(A.transpose(), Sw_At)))
            e = (b - np.matmul(self.weights, A)).transpose()
            self.weights = self.weights + np.matmul(K, e).transpose()

        err_data = np.matmul(self.weights, H) - pos_data
        train_err = np.linalg.norm(err_data, axis=1)

        return train_err

    def reconfig(self, n_kernels=None, kernels_std_scaling=None,
                 n_points=200, train_method='LS', end_points_constraints=False):

        if not n_kernels:
            n_kernels = self.n_kernels
        if not kernels_std_scaling:
            kernels_std_scaling = self.kernels_std_scaling

        s_data = np.linspace(0, 1, n_points)
        pos_data = np.hstack([self.get_pos(s) for s in s_data])
        # reconfigure MP
        self.__init_helper__(n_dofs=self.n_dofs, n_kernels=n_kernels, kernels_std_scaling=kernels_std_scaling)
        return self.train(s_data, pos_data, train_method=train_method, end_points_constraints=end_points_constraints)

    def to_state_dict(self):
        return {'weights': self.weights, 'config': self._config}

    @classmethod
    def from_state_dict(cls, state_dict):
        mp = cls(**state_dict['config'])
        mp.weights = state_dict['weights']
        return mp

    def save(self, filename):
        pickle.dump(self.to_state_dict(), open(filename, 'wb'))

    @staticmethod
    def load(filename):
        return MovementPrimitive.from_state_dict(pickle.load(open(filename, 'rb')))

    def deep_copy(self):
        return copy.deepcopy(self)

    def regress_vec(self, s: float) -> np.array:

        # take appropriate actions when x causes phi = 0 due to finite
        # numerical precision.
        if s < self.s_min:
            psi = np.zeros((self.n_kernels, 1))
            psi[0] = 1.0
        elif s > self.s_max:
            psi = np.zeros((self.n_kernels, 1))
            psi[-1] = 1.0
        else:
            psi = self._kernel_fun(s)

        phi = psi / np.sum(psi)
        return phi

    def regress_vec_dot(self, s: float, s_dot: float) -> np.array:

        if s < self.s_min or s > self.s_max:
            return np.zeros((self.n_kernels, 1))

        psi = self._kernel_fun(s)
        psi_dot = self._kernel_fun_dot(s, s_dot)
        sum_psi = np.sum(psi)
        sum_psi_dot = np.sum(psi_dot)

        phi = psi / sum_psi
        phi_dot = (psi_dot - phi * sum_psi_dot) / sum_psi
        return phi_dot

    def regress_vec_ddot(self, s: float, s_dot: float, s_ddot: float) -> np.array:

        if s < self.s_min or s > self.s_max:
            return np.zeros((self.n_kernels, 1))

        psi = self._kernel_fun(s)
        psi_dot = self._kernel_fun_dot(s, s_dot)
        psi_ddot = self._kernel_fun_ddot(s, s_dot, s_ddot)
        sum_psi = np.sum(psi)
        sum_psi_dot = np.sum(psi_dot)
        sum_psi_ddot = np.sum(psi_ddot)

        phi = psi / sum_psi
        phi_dot = (psi_dot - phi * sum_psi_dot) / sum_psi
        phi_ddot = (psi_ddot - 2 * phi_dot * sum_psi_dot - phi * sum_psi_ddot) / sum_psi
        return phi_ddot

    def _kernel_fun(self, s: float) -> np.array:
        return np.exp(-self.h * np.power(s - self.c, 2))

    def _kernel_fun_dot(self, s: float, s_dot: float) -> np.array:
        psi = self._kernel_fun(s)
        a = (s - self.c) * s_dot
        psi_dot = -2 * self.h * (psi * a)
        return psi_dot

    def _kernel_fun_ddot(self, s: float, s_dot: float, s_ddot: float) -> np.array:

        psi = self._kernel_fun(s)
        psi_dot = self._kernel_fun_dot(s, s_dot)
        a = (s - self.c) * s_dot
        a_dot = (s - self.c) * s_ddot + s_dot ** 2
        psi_ddot = -2 * self.h * (psi_dot * a + psi * a_dot)

        return psi_ddot


def mp_DPVP_training(path: np.array, vel_profile: np.array, n_kernels: Union[np.array, int]):
    if np.isscalar(n_kernels):
        n_kernels = [30, 200, n_kernels]

    if len(n_kernels) != 3:
        raise AttributeError('n_kernels must have length 3! Given input has length=%d' % len(n_kernels))

    n_dofs = path.shape[0]

    mp = MovementPrimitive(n_dofs=n_dofs, n_kernels=n_kernels[0])  # encodes the initial path
    mp2 = MovementPrimitive(n_dofs=n_dofs, n_kernels=n_kernels[1])  # encodes the up-sampled path
    mp3 = MovementPrimitive(n_dofs=n_dofs, n_kernels=n_kernels[2])  # encodes the path with the desired velocity profile

    # ========= Phase 1: train with traj-length, const vel =========
    sd_data = np.cumsum(np.concatenate([[0], np.linalg.norm(np.diff(path, axis=1), axis=0)]))
    sd_data = sd_data / sd_data[-1]  # normalize to be in [0, 1]

    ind = np.concatenate([[0], np.squeeze(np.where(np.diff(sd_data) > 1e-12)) + 1])
    sd0_data = sd_data[ind]
    Pd0_data = path[:, ind]

    Pd1_data = Pd0_data  # np.stack([np.interp(np.linspace(0, 1, 1000), sd0_data, Pd0_data[i]) for i in range(n_dofs)], axis=0)
    sd1_data = np.cumsum(np.concatenate([[0], np.linalg.norm(np.diff(Pd1_data, axis=1), axis=0)]))
    sd1_data = sd1_data / sd1_data[-1]
    train_err_1 = mp.train(sd1_data, Pd1_data, 'LS', end_points_constraints=False)
    # print('Phase 1 train_err:', train_err_1)

    # ========= Phase 2: Up-sample and re-train traj-length, constant vel =========
    Pd2_data = np.hstack([mp.get_pos(s) for s in np.linspace(0., 1., 1000)])  # up-sample
    s2_data = np.cumsum(np.concatenate([[0], np.linalg.norm(np.diff(Pd2_data, axis=1), axis=0)]))
    s2_data = s2_data / s2_data[-1]  # normalize to be in [0, 1]

    train_err_2 = mp2.train(s2_data, Pd2_data, 'LS', end_points_constraints=False)
    # print('Phase 2 train_err:', train_err_2)

    # ========= Phase 3-a: create new train data with desired velocity profile =========
    sd3_dot_data = vel_profile
    Time = np.linspace(0, 1, len(sd3_dot_data))
    sd3_data = np.concatenate([[0], np.cumsum(sd3_dot_data[:-1]*np.diff(Time))])
    vel_profile = vel_profile / sd3_data[-1]
    sd3_data = sd3_data / sd3_data[-1]
    Pd3_data = np.hstack([mp2.get_pos(s) for s in sd3_data])

    # ========= Phase 3-b: Train the MP to have the demo path and the desired velocity profile =========
    train_err_3 = mp3.train(Time, Pd3_data, 'LS', end_points_constraints=True)
    # print('Phase 3 train_err:', train_err_3)

    path_hat = np.hstack([mp3.get_pos(s) for s in np.linspace(0, 1, path.shape[1])])
    l = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(path_hat, axis=1), axis=0))])
    L = l[-1]
    l_q = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(path, axis=1), axis=0))])
    path_hat = np.stack([np.interp(l_q, l, path_hat[i]) for i in range(n_dofs)], axis=0)
    path_err = np.sum(np.linalg.norm(path_hat - path, axis=0)) / path.shape[1]

    # notice the use of s_dot = 1 / L  to produce a normalized velocity profile, whose integral equals 1.
    vel_profile_hat = np.linalg.norm(np.hstack([mp3.get_vel(s, 1.0 / L) for s in np.linspace(0, 1, len(vel_profile))]),
                                     axis=0)
    vel_err = np.linalg.norm(vel_profile_hat - vel_profile) / vel_profile.size

    # print('path_err = %f' % path_err)
    # print('vel_err = %f' % vel_err)

    # mp_vel = np.linalg.norm(np.hstack([mp.get_vel(s, 1.0 / L) for s in np.linspace(0, 1, len(vel_profile))]), axis=0)
    # mp2_vel = np.linalg.norm(np.hstack([mp2.get_vel(s, 1.0 / L) for s in np.linspace(0, 1, len(vel_profile))]), axis=0)

    # fig, ax = plt.subplots(1, 2)
    # ax[0].plot(path_hat[0], path_hat[1], color='magenta', linestyle='-', label='MP')
    # ax[0].plot(path[0], path[1], color='blue', linestyle='--', label='demo')
    # ax[0].set_title('Path')
    # # ax[1].plot(np.linspace(0, 1, len(mp_vel)), mp_vel, color=[0.85, 0.33, 0.1], linestyle='-', label='MP')
    # # ax[1].plot(np.linspace(0, 1, len(mp2_vel)), mp2_vel, color='cyan', linestyle='--', label='MP2')
    # ax[1].plot(np.linspace(0, 1, len(vel_profile_hat)), vel_profile_hat, color='magenta', linestyle='-', label='MP3')
    # ax[1].plot(np.linspace(0, 1, len(vel_profile)), vel_profile, color='blue', linestyle='--', label='demo')
    # ax[1].legend()
    # ax[1].set_title('Velocity profile')
    # plt.show()
    # input()
    # exit()

    return mp3, (path_err, vel_err)


def fifth_order_traj(t: float,
                     p0: Union[float, np.array],
                     pf: Union[float, np.array],
                     total_time: float) -> (Union[float, np.array], Union[float, np.array], Union[float, np.array]):

    if np.isscalar(p0):
        pos = p0
        vel = 0
        accel = 0
    else:
        n_dof = len(p0)
        pos = p0
        vel = np.zeros(n_dof)
        accel = np.zeros(n_dof)

    if t < 0:
        pos = p0
    elif t > total_time:
        pos = pf
    else:
        pos = p0 + (pf - p0) * (10 * pow(t / total_time, 3) -
                                15 * pow(t / total_time, 4) + 6 * pow(t / total_time, 5))
        vel = (pf - p0) * (30 * pow(t, 2) / pow(total_time, 3) -
                           60 * pow(t, 3) / pow(total_time, 4) + 30 * pow(t, 4) / pow(total_time, 5))
        accel = (pf - p0) * (60 * t / pow(total_time, 3) -
                             180 * pow(t, 2) / pow(total_time, 4) + 120 * pow(t, 3) / pow(total_time, 5))

    return pos, vel, accel


class DrawOnImage:

    def __init__(self, image_height=480, image_width=640, draw_type='trajectory',
                 record_stops=True, logger_Ts=0.02, dist_thres=-1, point_radius=5):

        self.win_name = 'Draw Path'

        cv2.namedWindow(winname=self.win_name)

        self.dist_thres = dist_thres  # register a new point if the  L2 distance in pixels from the last point is greater than this threshold
        self.point_radius = point_radius
        self.record_stops = record_stops
        self.logger_Ts = logger_Ts

        # trajectory params
        self.Time = []  # timestamp of each recorded point
        self.points = []

        self.color = (0, 0, 255)
        self.drawing = False
        self.img0 = np.zeros((image_height, image_width, 3), np.uint8)
        self.img = self.img0.copy()

        self.t0 = time.time()

        self.new_point = []
        self.lock = Lock()
        self.run_logger = True

        self.draw_type = draw_type.lower()

        if self.draw_type == 'line':
            self.points = [(-1, -1), (-1, -1)]
            cv2.setMouseCallback(self.win_name, DrawOnImage.__draw_line_callback, self)
        else:
            cv2.setMouseCallback(self.win_name, DrawOnImage.__draw_trajectory_callback, self)

    def __del__(self):
        cv2.destroyWindow(self.win_name)

    def __call__(self, img: np.array) -> Tuple[np.array, np.array, np.array]:

        if self.record_stops:
            self.run_logger = True
            logger_thread = Thread(target=self.logger_thread_fun)
            logger_thread.start()

        self.img0 = img.copy()
        self.__reset_image()
        while True:
            cv2.imshow(self.win_name, self.img)
            key = cv2.waitKey(1)

            if key == 27:  # Esc (focus must be on cv window)
                self.__reset_image()

            elif key == 13:  # 'enter'
                self.__reset_image()
                break

            elif key in (ord('u'), ord('U')):  # reset image
                self.__reset_image()

            elif key in (ord('l'), ord('L')):
                self.draw_type = 'line'
                cv2.setMouseCallback(self.win_name, DrawOnImage.__draw_line_callback, self)

            elif key in (ord('t'), ord('T')):
                self.draw_type = 'trajectory'
                cv2.setMouseCallback(self.win_name, DrawOnImage.__draw_trajectory_callback, self)

            elif key in (ord('p'), ord('P')):
                print('Draw type:', self.draw_type)

        if self.record_stops:
            self.run_logger = False
            logger_thread.join()

        x_data, y_data = map(np.array, zip(*self.points))
        # y_data = [-y for y in y_data]

        return x_data, y_data, np.array(self.Time)

    def __reset_image(self):
        self.img = self.img0.copy()

    def logger_thread_fun(self):

        while self.run_logger:
            if self.drawing:
                self.lock.acquire()
                self.points.append(self.new_point)
                self.Time.append(time.time() - self.t0)
                self.lock.release()
            time.sleep(self.logger_Ts)

    @staticmethod
    def __draw_line_callback(event, x, y, flags, self):

        if self.drawing:
            self.points[1] = (x, y)

        if event == cv2.EVENT_LBUTTONDOWN:
            self.drawing = True
            self.points[0] = (x, y)

        elif event == cv2.EVENT_MOUSEMOVE:
            if self.drawing:
                self.__reset_image()
                cv2.line(self.img, pt1=self.points[0], pt2=self.points[1], color=self.color, thickness=4)

        elif event == cv2.EVENT_LBUTTONUP:
            self.drawing = False
            self.__reset_image()
            cv2.line(self.img, pt1=self.points[0], pt2=self.points[1], color=self.color, thickness=4)

    @staticmethod
    def __draw_trajectory_callback(event, x, y, flags, self):

        if event == cv2.EVENT_LBUTTONDOWN:
            self.__reset_image()
            self.drawing = True

            self.lock.acquire()
            new_point = (x, y)
            self.points = [new_point]
            self.t0 = time.time()
            self.Time = [0.]
            cv2.circle(self.img, new_point, self.point_radius, self.color, thickness=-1)
            self.new_point = new_point
            self.lock.release()

        elif event == cv2.EVENT_LBUTTONUP:  # must be before the next elif, otherwise it will never stop!
            self.drawing = False

        elif event == cv2.EVENT_MOUSEMOVE and self.drawing:
            last_point = self.points[-1]
            new_point = (x, y)
            if np.linalg.norm(np.array(last_point) - np.array(new_point)) > self.dist_thres:
                self.lock.acquire()
                self.new_point = new_point
                self.lock.release()
                if not self.record_stops:
                    self.points.append(new_point)
                    self.Time.append(time.time() - self.t0)
                cv2.circle(self.img, new_point, self.point_radius, self.color, thickness=-1)


IM_HEIGHT = 500
IM_WIDTH = 500


def record_velocity_profile(save_filename='velocity_profile.pkl'):

    print('\33[1m\33[34m', 'Velocity Profile recording... (press [enter] to finish)', '\33[0m')

    xd_data, yd_data, Timed = hd_p(img)

    xd_data, yd_data = map(np.array, (xd_data, yd_data))
    yd_data = IM_HEIGHT - yd_data
    yd_data = yd_data - np.max([yd_data[0], yd_data[-1]])
    ind = np.argwhere(yd_data >= 0.).squeeze()
    yd_data = yd_data[ind]
    xd_data = xd_data[ind]
    yd_data = yd_data / np.max(yd_data)
    xd_data = xd_data - np.min(xd_data)
    sd_data = xd_data / np.max(xd_data)

    mp = MovementPrimitive(n_dofs=1, n_kernels=min(100, len(sd_data)//2))
    mp.train(sd_data, yd_data, train_method='LS', end_points_constraints=True)

    s_data = np.linspace(0, 1, 1000)
    y_data = np.hstack([mp.get_pos(s) for s in s_data])

    velocity_profile = np.squeeze(y_data)

    if save_filename:
        pickle.dump(velocity_profile, open(save_filename, 'wb'))

    # velocity_profile = pickle.load(open(filename, 'rb'))
    # fig, ax = plt.subplots()
    # ax.plot(sd_data, yd_data.squeeze(), label='demo')
    # ax.plot(s_data, velocity_profile.squeeze(), label='mp')
    # ax.legend()
    # plt.show()
    # input()

    return velocity_profile


def record_demo(save_filename=None):

    print('\33[1m\33[34m', 'Path recording... (press [enter] to finish)', '\33[0m')

    xd_data, yd_data, Timed = hd_p(img)

    yd_data = IM_HEIGHT - yd_data

    xd_data = xd_data / IM_WIDTH
    yd_data = yd_data / IM_HEIGHT

    xd_data = xd_data - xd_data[0]
    yd_data = yd_data - yd_data[0]

    if save_filename:
        pickle.dump((xd_data, yd_data, Timed), open(save_filename, 'wb'))

    return xd_data, yd_data, Timed


if __name__ == '__main__':

    data_filename = sys.argv[1]

    plt.ion()

    hd_p = DrawOnImage(IM_HEIGHT, IM_WIDTH, draw_type='trajectory',
                       record_stops=True, logger_Ts=0.02, dist_thres=-1, point_radius=4)

    img = np.zeros((IM_HEIGHT, IM_WIDTH, 3), np.uint8)

    xd_data, yd_data, Timed = record_demo()
    demo_data = np.stack([xd_data, yd_data], axis=0)

    vel_profile = record_velocity_profile()

    fig, ax = plt.subplots()
    ax.plot(vel_profile)
    ax.set_ylabel('recorded velocity profile')

    # vel_profile = np.array([fifth_order_traj(t, 0.0, 1.0, 1.0)[1] for t in np.linspace(0, 1, 1000)])

    mp, _ = mp_DPVP_training(path=demo_data, vel_profile=vel_profile, n_kernels=40)

    s_data = np.linspace(0, 1, 200)
    pos = np.hstack([mp.get_pos(s) for s in s_data])
    vel = np.hstack([mp.get_vel(s, 1.0) for s in s_data])
    vel_norm = np.linalg.norm(vel, axis=0)

    fig, ax = plt.subplots(2, 1)
    ax[0].plot(pos[0, :], pos[1, :])
    ax[0].set_ylabel('2D path')
    ax[1].plot(s_data, vel_norm)
    ax[1].set_ylabel('vel norm')
    plt.pause(0.001)

    scipy.io.savemat(data_filename, {'sd_data': s_data, 'Pd_data': pos})

    input('Press [enter] to exit...')
