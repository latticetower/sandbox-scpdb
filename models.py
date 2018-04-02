""

from keras.models import Sequential
from keras.layers import Conv3D, MaxPooling3D, Dropout, Flatten, Dense, Activation
from keras.layers.advanced_activations import ELU
from keras.optimizers import Adam


def get_deepsite_model(boxsize=16):
    model = Sequential()
    model.add(Conv3D(32, (8, 8, 8), padding='same', input_shape=(boxsize, boxsize, boxsize, 8),
                    data_format='channels_last'))
    model.add(ELU())
    model.add(Conv3D(48, (4, 4, 4)))
    model.add(ELU())
    model.add(MaxPooling3D(pool_size=(2, 2, 2)))
    model.add(Dropout(0.25))
    model.add(Conv3D(64, (4, 4, 4), padding='same'))
    model.add(ELU())
    model.add(Conv3D(96, (4, 4, 4)))
    model.add(ELU())
    model.add(MaxPooling3D(pool_size=(2, 2, 2)))
    model.add(Dropout(0.25))
    model.add(Flatten())
    model.add(Dense(256))
    model.add(ELU())
    model.add(Dropout(0.5))
    model.add(Dense(1))
    model.add(Activation('sigmoid'))
    opt = Adam()
    model.compile(loss='binary_crossentropy', optimizer=opt, metrics=['accuracy'])
    return model